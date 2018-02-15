module assimilation.likelihood.DiscreteExperimentalLikelihood;

import std.algorithm;
import std.array;
import std.conv;
import std.datetime;
import std.math;
import std.parallelism;
import std.range;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import data.Ensemble;
import data.Timeseries;
import data.Vector;
import integrators.Integrator;
import utility.Normal;

class DiscreteExperimentalLikelihood : LikelihoodGetter {

    Integrator integrator;
    double minimumOffset; ///The most a true time can be less than the errant time; this is equal to the most an errant time can be more than the truth
    double maximumOffset; ///The most a true time can be more than the errant time; this is equal to the most an errant time can be less than the truth
    uint bins; ///The amount of bins into which to sort the time likelihood
    int[] timeLikelihood;

    /**
     * Constructs a likelihood getter with information about the experiment, as well as a priori knowledge of time offset and desired number of bins for time offset likelihood
     * Also integrates ensemble timeseries to the maximum offset in order to have values for the interval from minimum to maximum offset
     */
    this(Timeseries!Vector observations, Vector stateError, Integrator integrator, double minimumOffset, double maximumOffset, uint bins) {
        super(observations, stateError);
        this.integrator = integrator;
        this.minimumOffset = minimumOffset;
        this.maximumOffset = maximumOffset;
        this.bins = bins;
        foreach(i; 0..bins) {
            this.timeLikelihood ~= 0;
        }
    }

    /**
     * Returns likelihood packaged with discretely defined experimentally determined for a given time
     */
    override Likelihood opCall(double time, Timeseries!Ensemble ensembles) {
        Ensemble ensemble = ensembles.members[$ - 1];
        assert(ensembles.times[$ - 1] + maximumOffset > ensembles.times[$ - 1]);
        //The following code does not work because of dlang iota implementation, so I have bypassed that manually for now
        /*foreach(i; iota(start, end, step)) {
            ensemble = this.integrator(ensemble, ensembles.dt);
            ensembles.add(i + ensembles.dt, ensemble);
        }*/
        double start = ensembles.times[$ - 1];
        double end = start + this.maximumOffset;
        double step = ensembles.dt;
        int steps = cast(int)((end - start) / step);
        foreach(i; 0..steps) {
            double placeTime = start + step * i;
            ensemble = this.integrator(ensemble, step);
            ensembles.add(placeTime + step, ensemble);
        }
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        ensemble = ensembles.value(time, this.integrator);
        assert(ensembles !is null);
        int[] timeLikelihood = this.getTimeLikelihood(time, ensembles);
        Vector[] observationTrajectory = this.getObservationTrajectory(time);
        assert(timeLikelihood.length == observationTrajectory.length, "timeLikelihood has " ~ timeLikelihood.length.to!string ~ " elements whereas observationTrajectory has " ~ observationTrajectory.length.to!string);
        //Apply likelihoods
        double[] xLikelihood;
        double[] yLikelihood;
        double[] zLikelihood;
        foreach(i; 0..ensemble.size) {
            xLikelihood ~= 0;
            yLikelihood ~= 0;
            zLikelihood ~= 0;
        }
        foreach(index, ref component; observationTrajectory.parallel) {
            xLikelihood = xLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.x, this.stateError.x)).array;
            yLikelihood = yLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.y, this.stateError.y)).array;
            zLikelihood = zLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.z, this.stateError.z)).array;
        }
        double xSum = xLikelihood.sum;
        double ySum = yLikelihood.sum;
        double zSum = zLikelihood.sum;
        //writeln("Raw Sums ", xSum, " ", ySum, " ", zSum);
        //If there is an issue with sums normalize
        //This code is unsustainable so don't use it to solve the problem;
        if(xSum == 0 || ySum == 0 || zSum == 0) {
            xLikelihood = xLikelihood.map!(a => a + 0.1).array;
            yLikelihood = yLikelihood.map!(a => a + 0.1).array;
            zLikelihood = zLikelihood.map!(a => a + 0.1).array;
        }
        xSum = xLikelihood.sum;
        ySum = yLikelihood.sum;
        zSum = zLikelihood.sum;
        foreach(index, ref component; timeLikelihood.parallel) {
            xLikelihood[index] /= xSum;
            yLikelihood[index] /= ySum;
            zLikelihood[index] /= zSum;
        }
        //writeln("Likelihood Sums: ", xLikelihood.sum, " ", yLikelihood.sum, " ", zLikelihood.sum);
        return new Likelihood(xLikelihood, yLikelihood, zLikelihood);
    }

    /**
     * At a certain time, determines discretely defined likelihood in time
     */
    int[] getTimeLikelihood(double time, Timeseries!Ensemble ensembles) {
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        Vector obs = this.observations.value(time);
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        double[] binMiddles = iota(time + this.minimumOffset + binWidth / 2, time + this.maximumOffset, binWidth).array;
        //import std.stdio;
        //writeln("bin times ", binMiddles);
        //writeln("Ensemble mean ", ensembles.members[$-1].eMean);
        Vector[] binValues = binMiddles.map!(a => ensembles.meanSeries.value(a, this.integrator)).array;
        int[] binQuantities; 
        foreach(i; 0..bins) {
            binQuantities ~= 0;
        }
        Vector distance;
        foreach(i; 0..binMiddles.length) {
            distance = binValues[i] - obs;
            //writeln("Distance: ", binValues[i], " - ", obs, " = ", distance);
            //writeln("State error: ", this.stateError);
            if(abs(distance.x) < 2 * this.stateError.x && abs(distance.y) < 2 * this.stateError.y && abs(distance.z) < 2 * this.stateError.z) {
                binQuantities[i] += 1;
            }
        }
        foreach(index, ref component; binQuantities.parallel) {
            this.timeLikelihood[index] += component;
        }
        return this.timeLikelihood;
    }

    /**
     * Performs getTimeLikelihood, but returns time likelihood normalized to a discrete PDF
     */
    double[] getNormalizedTimeLikelihood(double time, Timeseries!Ensemble ensembles) {
        double[] normalizedLikelihood = cast(double[])this.getTimeLikelihood(time, ensembles);
        immutable double sum = normalizedLikelihood.sum;
        foreach(ref component; normalizedLikelihood.parallel) {
            component /= sum;
        }
        return normalizedLikelihood;
    }

    /**
     * Gets baselines for each bin to find likelihood for a given observation 
     * TODO: don't integrate backwards
     */
    Vector[] getObservationTrajectory(double time) {
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        Vector obs = this.observations.value(time);
        Vector[] baselines;
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        Vector binEdge = this.integrator.integrateTo(obs, this.minimumOffset, 10 * cast(uint)abs(this.minimumOffset) / this.bins);
        Vector state = this.integrator.integrate(binEdge, binWidth / 2);
        baselines ~= state;
        foreach(i; 0..this.bins - 1) {
            state = this.integrator.integrate(binEdge, binWidth);
            baselines ~= state;
        }
        return baselines;
    }

}