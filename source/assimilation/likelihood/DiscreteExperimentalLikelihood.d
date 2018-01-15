module assimilation.likelihood.DiscreteExperimentalLikelihood;

import std.algorithm;
import std.array;
import std.math;
import std.parallelism;
import std.range;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import data.Ensemble;
import data.Timeseries;
import data.Vector;
import integrators.Integrator;

class DiscreteExperimentalLikelihood : LikelihoodGetter {

    Timeseries!Ensemble ensembles;
    Integrator integrator;
    double minimumOffset; ///The most a true time can be less than the errant time; this is equal to the most an errant time can be more than the truth
    double maximumOffset; ///The most a true time can be more than the errant time; this is equal to the most an errant time can be less than the truth
    ulong bins; ///The amount of bins into which to sort the time likelihood
    int[] timeLikelihood;

    /**
     * Constructs a likelihood getter with information about the experiment, as well as a priori knowledge of time offset and desired number of bins for time offset likelihood
     * Also integrates ensemble timeseries to the maximum offset in order to have values for the interval from minimum to maximum offset
     */
    this(Timeseries!Vector observations, Timeseries!Ensemble ensembles, Vector stateError, Integrator integrator, double minimumOffset, double maximumOffset, ulong bins) {
        super(observations, stateError);
        this.ensembles = ensembles;
        this.integrator = integrator;
        this.minimumOffset = minimumOffset;
        this.maximumOffset = maximumOffset;
        this.bins = bins;
        this.timeLikelihood[0..bins] = 0;
        Ensemble ensemble = this.ensembles.members[$];
        foreach(i; iota(this.ensembles.times[$], this.ensembles.times[$] + maximumOffset, this.ensembles.dt)) {
            ensemble = this.integrator(ensemble, this.ensembles.dt);
            this.ensembles.add(i + this.ensembles.dt, ensemble);
        }
    }

    /**
     * Returns likelihood packaged with discretely defined Gaussian likelihoods for a given time
     */
    override Likelihood opCall(double time) {
        assert(this.observations.times.canFind(time));
        Ensemble ensemble = ensembles.value(time, this.integrator);
        return null;
    }

    /**
     * At a certain time, determines discretely defined likelihood in time
     */
    int[] getTimeLikelihood(double time) {
        assert(this.observations.times.canFind(time));
        Vector obs = this.observations.value(time);
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        double[] binMiddles = iota(this.minimumOffset + binWidth / 2, this.maximumOffset, binWidth).array;
        Vector[] binValues = binMiddles.map!(a => this.ensembles.meanSeries.value(a)).array;
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
    double[] getNormalizedTimeLikelihood(double time) {
        double[] normalizedLikelihood = cast(double[])this.getTimeLikelihood(time);
        immutable double sum = normalizedLikelihood.sum;
        foreach(ref component; normalizedLikelihood.parallel) {
            component /= sum;
        }
        return normalizedLikelihood;
    }

}