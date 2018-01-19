module assimilation.likelihood.DiscreteExperimentalLikelihood;

import std.algorithm;
import std.array;
import std.conv;
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
        Ensemble ensemble = ensembles.value(ensembles.length - 1);
        foreach(i; iota(ensembles.times[$ - 1], ensembles.times[$ - 1] + maximumOffset, ensembles.dt)) {
            ensemble = this.integrator(ensemble, ensembles.dt);
            ensembles.add(i + ensembles.dt, ensemble);
        }
        assert(this.observations.times.canFind(time));
        ensemble = ensembles.value(time, this.integrator);
        int[] timeLikelihood = this.getTimeLikelihood(time, ensembles);
        Vector[] observationTrajectory = this.getObservationTrajectory(time);
        assert(timeLikelihood.length == observationTrajectory.length, "timeLikelihood has " ~ timeLikelihood.length.to!string ~ " elements whereas observationTrajectory has " ~ observationTrajectory.length.to!string);
        //Apply likelihoods
        double[] xLikelihood;
        xLikelihood[0..cast(uint)ensemble.size] = 0.0;
        double[] yLikelihood;
        yLikelihood[0..cast(uint)ensemble.size] = 0.0;
        double[] zLikelihood;
        zLikelihood[0..cast(uint)ensemble.size] = 0.0;
        foreach(index, ref component; observationTrajectory.parallel) {
            xLikelihood = xLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.x, this.stateError.x)).array;
            yLikelihood = yLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.y, this.stateError.y)).array;
            zLikelihood = zLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.z, this.stateError.z)).array;
        }
        immutable double xSum = xLikelihood.sum;
        immutable double ySum = yLikelihood.sum;
        immutable double zSum = zLikelihood.sum;
        foreach(index, ref component; timeLikelihood.parallel) {
            xLikelihood[index] /= xSum;
            yLikelihood[index] /= ySum;
            zLikelihood[index] /= zSum;
        }
        return new Likelihood(xLikelihood, yLikelihood, zLikelihood);
    }

    /**
     * At a certain time, determines discretely defined likelihood in time
     */
    int[] getTimeLikelihood(double time, Timeseries!Ensemble ensembles) {
        assert(this.observations.times.canFind(time));
        Vector obs = this.observations.value(time);
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        double[] binMiddles = iota(this.minimumOffset + binWidth / 2, this.maximumOffset, binWidth).array;
        Vector[] binValues = binMiddles.map!(a => ensembles.meanSeries.value(a)).array;
        int[] binQuantities; 
        binQuantities[0..this.bins] = 0;
        Vector distance;
        foreach(i; 0..binMiddles.length) {
            distance = binValues[i] - obs;
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
        assert(this.observations.times.canFind(time));
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