module assimilation.likelihood.DiscreteExperimentalLikelihood;

import std.algorithm;
import std.array;
import std.conv;
import std.datetime;
import std.math;
import std.parallelism;
import std.range;
import mir.random;
import mir.random.variable;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import data.Ensemble;
import data.Timeseries;
import data.Vector;
import integrators.Integrator;
import math.Matrix;
import utility.ArrayStats;
import utility.Normal;

/**
 * Gets likelihood inferentially using Bayes' rule
 */
class DiscreteExperimentalLikelihood : LikelihoodGetter {

    Integrator integrator;
    double minimumOffset; ///The most a true time can be less than the errant time; this is equal to the most an errant time can be more than the truth
    double maximumOffset; ///The most a true time can be more than the errant time; this is equal to the most an errant time can be less than the truth
    uint bins; ///The amount of bins into which to sort the time likelihood
    double[] timeLikelihood; ///Histogram of time likelihood
    double timeOffset; ///The known offset in time; only use if testing knownErrorNormalLikelihood
    double timeError; ///The known standard deviation of time error; ^^^

    /**
     * Gets the mathematical expected value for time offset. This is the weighted average of the bin middles.
     */
    @property double expectedTime() {
        assert(bins == timeLikelihood.length, "Likelihood length is not the same as bin length.");
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        double[] binMiddles = [this.minimumOffset + binWidth / 2];
        foreach(i; 1..bins) {
            binMiddles ~= binMiddles[i - 1] + binWidth;
        }
        double[] likelihood = this.normalizedTimeLikelihood;
        foreach(index, ref component; binMiddles.parallel) {
            component *= likelihood[index];
        } 
        return binMiddles.sum;
    }

    /**
     * Returns standard deviation of time likelihood fit to a Gaussian
     */
    @property double timeDeviation() {
        double expectedOffset = this.expectedTime;
        assert(bins == timeLikelihood.length);
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        double[] binMiddles = [this.minimumOffset + binWidth / 2];
        foreach(i; 1..bins) {
            binMiddles ~= binMiddles[i - 1] + binWidth;
        }  
        double[] likelihood = this.normalizedTimeLikelihood;
        foreach(index, ref component; binMiddles.parallel) {
            component = ((component - expectedOffset).pow(2)) * likelihood[index];
        } 
        return sqrt(binMiddles.sum);
    }

    /** 
     * Performs a chi-square test on the discrete likelihood distribution with its normal fit
     */
    @property double timeGaussianity() {
        double expectedMean = this.expectedTime;
        double expectedDeviation = this.timeDeviation;
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        double[] binMiddles = [this.minimumOffset + binWidth / 2];
        double[] likelihood = this.normalizedTimeLikelihood;
        double[] expectedValues = binMiddles.dup.map!(a => normalVal(a, expectedMean, expectedDeviation)).array;
        double expectedSum = expectedValues.sum;
        foreach(ref component; expectedValues) {
            component /= expectedSum;
        }
        foreach(index, ref component; likelihood) {
            component = ((component - expectedValues[index]).pow(2)) / expectedValues[index];
        }
        return likelihood.sum;
    }

    /**
     * Returns time likelihood normalized to a discrete PDF
     */
    @property double[] normalizedTimeLikelihood() {
        double[] normalizedLikelihood = this.timeLikelihood.to!(double[]).dup;
        immutable double sum = normalizedLikelihood.sum;
        foreach(ref component; normalizedLikelihood.parallel) {
            component /= sum;
        }
        return normalizedLikelihood;
    }

    /**
     * Constructs a likelihood getter with information about the experiment, as well as a priori knowledge of time offset and desired number of bins for time offset likelihood
     * Also integrates ensemble timeseries to the maximum offset in order to have values for the interval from minimum to maximum offset
     */
    this(Timeseries!Vector observations, Vector stateError, Integrator integrator, double minimumOffset, double maximumOffset, uint bins, double timeOffset = 0, double timeError = 0) {
        super(observations, stateError);
        this.integrator = integrator;
        this.minimumOffset = minimumOffset;
        this.maximumOffset = maximumOffset;
        this.bins = bins;
        foreach(i; 0..bins) {
            this.timeLikelihood ~= 0;
        }
        this.timeLikelihood[bins / 2 - 1] = 1;
        this.timeLikelihood[bins / 2] = 1;
        this.timeOffset = timeOffset;
        this.timeError = timeError;
    }

    /**
     * Returns likelihood
     */
    override Likelihood opCall(double time, Timeseries!Ensemble ensembles) {
        //Choose one:
        //return this.likelihoodFromNormalTime(time, ensembles); //Only works with RHF for now
        //return this.likelihoodFromDiscreteTime(time, ensembles); //Only works with RHF for now
        return this.normalLikelihood(time, ensembles); //Only works with EAKF for now
        //return this.knownErrorNormalLikelihood(time, ensembles); //Only works with EAKF for now; requires minor changes to code for now
    }

    /** 
     * Fast alternative to normalLikelihood if error is known
     */
    Likelihood knownErrorNormalLikelihood(double time, Timeseries!Ensemble ensembles) {
        Vector obs = this.observations.value(time);
        double[] xPseudoMeasurements;
        double[] yPseudoMeasurements;
        double[] zPseudoMeasurements;            
        auto gen = Random(unpredictableSeed);
        foreach(i; 0..100) {
            auto newTime = NormalVariable!double(this.timeOffset, this.timeError);
            Vector base = this.integrator.integrateTo(obs, newTime(gen), 1);
            auto normalX = NormalVariable!double(base.x, this.stateError.x);
            auto normalY = NormalVariable!double(base.y, this.stateError.y);
            auto normalZ = NormalVariable!double(base.z, this.stateError.z);
            xPseudoMeasurements ~= normalX(gen);
            yPseudoMeasurements ~= normalY(gen);
            zPseudoMeasurements ~= normalZ(gen);
        }
        Vector meanVector = Vector(
            xPseudoMeasurements.mean,
            yPseudoMeasurements.mean,
            zPseudoMeasurements.mean
        );
        Vector deviationVector = Vector(
            standardDeviation!1(xPseudoMeasurements),
            standardDeviation!1(yPseudoMeasurements),
            standardDeviation!1(zPseudoMeasurements)  
        );
        return new Likelihood(meanVector, deviationVector);
    }

    /**
     * Returns a normal likelihood from normal-assumed time
     * TODO: Move into its own LikelihoodGetter
     */
    Likelihood normalLikelihood(double time, Timeseries!Ensemble ensembles) {
        Ensemble ensemble = new Ensemble(ensembles.members[$ - 1].members);
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
        this.getTimeLikelihood(time, ensembles);
        Vector obs = this.observations.value(time);
        double expectedOffset = this.expectedTime;
        double timeDeviation = this.timeDeviation;
        double[] xPseudoMeasurements;
        double[] yPseudoMeasurements;
        double[] zPseudoMeasurements;            
        auto gen = Random(unpredictableSeed);
        foreach(i; 0..100) {
            auto newTime = NormalVariable!double(expectedOffset, timeDeviation);
            Vector base = this.integrator.integrateTo(obs, newTime(gen), 1);
            auto normalX = NormalVariable!double(base.x, this.stateError.x);
            auto normalY = NormalVariable!double(base.y, this.stateError.y);
            auto normalZ = NormalVariable!double(base.z, this.stateError.z);
            xPseudoMeasurements ~= normalX(gen);
            yPseudoMeasurements ~= normalY(gen);
            zPseudoMeasurements ~= normalZ(gen);
        }
        Vector meanVector = Vector(
            xPseudoMeasurements.mean,
            yPseudoMeasurements.mean,
            zPseudoMeasurements.mean
        );
        Vector deviationVector = Vector(
            standardDeviation!1(xPseudoMeasurements),
            standardDeviation!1(yPseudoMeasurements),
            standardDeviation!1(zPseudoMeasurements)  
        );
        return new Likelihood(meanVector, deviationVector);
    }

    /**
     * Returns likelihood packaged with discretely defined experimentally determined pdf for a given time
     * Fits time likelihood to a Gaussian, then uses a Monte Carlo-style kernel density approximation with kernel widths 
     * given by manufacturer
     * Assumes time likelihood is Gaussian (data suggest it will be with current method; it is advised to use likelihoodFromDiscrete time otherwise)
     */
    Likelihood likelihoodFromNormalTime(double time, Timeseries!Ensemble ensembles) {
        Ensemble ensemble = new Ensemble(ensembles.members[$ - 1].members);
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
        this.getTimeLikelihood(time, ensembles);
        double[] xLikelihood;
        double[] yLikelihood;
        double[] zLikelihood;
        foreach(i; 0..ensemble.size) {
            xLikelihood ~= 0.0000001;
            yLikelihood ~= 0.0000001;
            zLikelihood ~= 0.0000001;
        }
        Vector obs = this.observations.value(time);
        double expectedOffset = this.expectedTime;
        double timeDeviation = this.timeDeviation;
        double[] pseudoTimes;
        //Generates pseudotimes; number of kernels can be changed within the code, but should probably be optimized
        auto gen = Random(unpredictableSeed);
        foreach(i; 0..100) {
            auto newTimeVar = NormalVariable!double(expectedOffset, timeDeviation);
            pseudoTimes ~= clamp(newTimeVar(gen), this.minimumOffset, this.maximumOffset);
        }
        Vector[] bases = pseudoTimes.dup.map!(a => Vector(a, 0, 0)).array;
        foreach(ref component; bases.parallel) {
            component = this.integrator.integrateTo(obs, component.x, 1);
            foreach(index, ref pointProbability; xLikelihood.parallel) {
                pointProbability += normalVal(ensemble.members[index].x, component.x, this.stateError.x);
                yLikelihood[index] += normalVal(ensemble.members[index].y, component.y, this.stateError.y);
                zLikelihood[index] += normalVal(ensemble.members[index].z, component.z, this.stateError.z);
            }
        }
        double xSum = xLikelihood.sum;
        double ySum = yLikelihood.sum;
        double zSum = zLikelihood.sum;
        assert(xSum != 0 && ySum != 0 && zSum != 0, "Experimental likelihood sum is 0");
        foreach(index, ref component; timeLikelihood.parallel) {
            xLikelihood[index] /= xSum;
            yLikelihood[index] /= ySum;
            zLikelihood[index] /= zSum;
        }
        return new Likelihood(xLikelihood, yLikelihood, zLikelihood);
    }

    /**
     * Returns likelihood packaged with discretely defined experimentally determined pdf for a given time
     * Assume likelihood histogram is accurate
     */
    Likelihood likelihoodFromDiscreteTime(double time, Timeseries!Ensemble ensembles) {
        Ensemble ensemble = new Ensemble(ensembles.members[$ - 1].members);
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
        double[] timeLikelihood = this.getTimeLikelihood(time, ensembles);
        Vector[] observationTrajectory = this.getObservationTrajectory(time);
        assert(timeLikelihood.length == observationTrajectory.length, "timeLikelihood has " ~ timeLikelihood.length.to!string ~ " elements whereas observationTrajectory has " ~ observationTrajectory.length.to!string);
        //Apply likelihoods
        double[] xLikelihood;
        double[] yLikelihood;
        double[] zLikelihood;
        foreach(i; 0..ensemble.size) {
            xLikelihood ~= 0.0000001;
            yLikelihood ~= 0.0000001;
            zLikelihood ~= 0.0000001;
        }
        foreach(index, ref component; observationTrajectory.parallel) {
            xLikelihood = xLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.x, this.stateError.x)).array;
            yLikelihood = yLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.y, this.stateError.y)).array;
            zLikelihood = zLikelihood.map!(a => a + timeLikelihood[index] * normalVal(a, component.z, this.stateError.z)).array;
        }
        double xSum = xLikelihood.sum;
        double ySum = yLikelihood.sum;
        double zSum = zLikelihood.sum;
        assert(xSum != 0 && ySum != 0 && zSum != 0, "Experimental likelihood sum is 0");
        foreach(index, ref component; timeLikelihood.parallel) {
            xLikelihood[index] /= xSum;
            yLikelihood[index] /= ySum;
            zLikelihood[index] /= zSum;
        }
        return new Likelihood(xLikelihood, yLikelihood, zLikelihood);
    }

    /**
     * Gets time likelihood, accounting for ensemble variance
     */
    double[] getTimeLikelihood(double time, Timeseries!Ensemble ensembles) {
        //Ensure that the observation time is close enough to one that we know
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        //Get observation at time
        Vector obs = this.observations.value(time);
        //Split the interval into bins; TODO: Perhaps make this a property
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        //Get the times in the centers of each bin
        double[] binMiddles = iota(time + this.minimumOffset + binWidth / 2, time + this.maximumOffset, binWidth).array;
        //Get the system state at the middle of each bin
        //Ensembles passed into function should already have been advanced to the end of the interval
        Vector[] binValues = binMiddles.map!(a => ensembles.meanSeries.value(a, this.integrator)).array;
        //Get the values of the ensemble points at the bin middles;
        //Might be faster to do this, then compute means from ensembles
        Ensemble[] ensembleValues = binMiddles.map!(a => ensembles.value(a, this.integrator)).array;
        //Initialize an array to store time likelihood - uniformly zero for now
        double[] binQuantities; 
        //Can maybe do this faster or more elegantly
        foreach(i; 0..this.bins) {
            binQuantities ~= 0;
        }
        foreach(i; 0..binMiddles.length) {
            Vector distance = binValues[i] - obs;
            //Be careful; options 1 and 2 are biased because they do not account for ensemble variance
            //However, they are computationally faster since option 3 involves matrix inversion and multiplication

            //Option 1: Discrete binary
            //If within 3 standard deviations, use it
            /*if(abs(distance.x) < 3 * this.stateError.x && abs(distance.y) < 3 * this.stateError.y && abs(distance.z) < 3 * this.stateError.z) {
                binQuantities[i] += 1;
            }*/

            //Option 2: Smooth
            //Assume error is dimensionally independent, then find probability by multiplying the probabilities
            /*binQuantities[i] = normalVal(distance.x, 0, this.stateError.x) 
                              * normalVal(distance.y, 0, this.stateError.y)
                              * normalVal(distance.z, 0, this.stateError.z);*/

            //Option 3: Smooth multivariate
            //Assume error is dimensionally independent, then find probability
            //By finding the value of a multivariate (in this case trivariate) normal pdf at the value of distance
            //The observation error covariance is 0 outside the diagonal because we assume obs error is dimensionally independent
            Matrix!(double, 3, 3) obsErrorCovariance = new Matrix!(double, 3, 3)([
                                                                                [this.stateError.x, 0, 0],
                                                                                [0, this.stateError.y, 0],
                                                                                [0, 0, this.stateError.z]
                                                                                ]);
            //The covariance matrix of the ensemble has the variance across the diagonals and the covariance everywhere else
            Matrix!(double, 3, 3) ensembleCovariance = covariance!(3, 1)([ensembleValues[i].xValues, ensembleValues[i].yValues, ensembleValues[i].zValues]);
            //The value of the probability is given by the trivariate normal PDF at the observation with mean ensemble mean
            //and covariance the sum of the ensemble covariance and the obs error covariance
            //This is the probability of taking the observation if the truth is taken from multivariate normal distribution
            //that is represented by the ensemble
            binQuantities[i] = trivariateNormalVal(obs.handle, binValues[i].handle, obsErrorCovariance + ensembleCovariance);
        }
        //In parallel, update the time likelihood with the new inferred likelihood
        foreach(index, ref component; binQuantities.parallel) {
            this.timeLikelihood[index] += component;
        }
        return this.timeLikelihood;
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
        Vector binEdge = this.integrator.integrateTo(obs, this.minimumOffset, 10 * cast(uint) (abs(this.minimumOffset) / binWidth));
        Vector state = this.integrator.integrate(binEdge, binWidth / 2);
        baselines ~= state;
        foreach(i; 0..this.bins - 1) {
            state = this.integrator.integrate(state, binWidth);
            baselines ~= state;
        }
        return baselines;
    }

}