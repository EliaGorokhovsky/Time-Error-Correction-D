/**
 * The source of the treatment for this experiment
 * A novel method of inflating likelihood by accounting for time error
 * This module's purpose is to estimate the magnitude of time error and then apply that
 * to likelihood, returning an updated one for assimilation
 */
module assimilation.likelihood.DiscreteExperimentalLikelihood;

import std.stdio;

import std.algorithm; //Used for map and reduce functions (list operations)
import std.array; //Used to convert map and reduce function outputs to arrays
import std.conv; //Used to cast things
import std.math; //Used to check approximate euality among doubles
import std.parallelism; //Used to run operations in parallel for speed
import std.range; //Used to generate iotas (ranges of numbers with steps)
import mir.random; //Used for random number generation
import mir.random.variable; //Used for random number generation
import assimilation.likelihood.Likelihood; //Used for outputs from this object's opCall
import assimilation.likelihood.LikelihoodGetter; //Used for the parent class of this object
import data.Ensemble; //Used to handle ensembles for inference
import data.Timeseries; //Used to handle ensembles and observations for inference
import integrators.Integrator; //Used to get values of timeseries when they are not defined
import math.Matrix; //Used to do covariance computations when finding probability that observation was taken at a given time
import math.Vector; //Used for data storage
import utility.ArrayStats; //Used for mean, standard deviation, and covariance
import utility.Normal; //Used to get the value of the normal probability density function

/**
 * Gets likelihood inferentially using Bayes' rule
 * Takes dimensionality of system being used
 */
class DiscreteExperimentalLikelihood(uint dim) : LikelihoodGetter!dim {

    Integrator!dim integrator;
    uint stepsTaken = 0; ///How many steps the algorithm has taken
    double minimumOffset; ///The most a true time can be less than the errant time; this is equal to the most an errant time can be more than the truth
    double maximumOffset; ///The most a true time can be more than the errant time; this is equal to the most an errant time can be less than the truth
    uint bins; ///The amount of bins into which to sort the time likelihood
    double timeOffset; ///The mean of the time error distribution - known a priori or inferred
    double timeError; ///The variance of the time error distribution - known a priori or inferred
    Random* gen; ///The random number generator used for the random generation; this is passed in so we can control the seed globally
    bool multiply = false; ///The inferential method can either multiply or add successive inferences; if true, multiply
    bool correctSlope = true; ///The inferential method can correct time error using knowledge about the system; if true, do so

    /**
     * Gets the standard deviation of the time error distribution
     */
    @property double timeDeviation() {
        return sqrt(this.timeError);
    }

    /**
     * Constructs a likelihood getter with information about the experiment, as well as a priori knowledge of time offset and desired number of bins for time offset likelihood
     * Also integrates ensemble timeseries to the maximum offset in order to have values for the interval from minimum to maximum offset
     */
    this(Timeseries!(Vector!(double, dim)) observations, Vector!(double, dim) stateError, Integrator!dim integrator, double minimumOffset, double maximumOffset, uint bins, Random* gen, double timeOffset = 0, double timeError = 0) {
        super(observations, stateError);
        this.integrator = integrator;
        this.minimumOffset = minimumOffset;
        this.maximumOffset = maximumOffset;
        this.bins = bins;
        this.timeOffset = timeOffset;
        this.timeError = timeError;
        this.gen = gen;
        //File("data/sandbox/test.csv", "a").writeln("Offset, Variance");
    }

    /**
     * Returns likelihood
     */
    override Likelihood!dim opCall(double time, Timeseries!(Ensemble!dim) ensembles) {
        //Choose one:
        //-------------------------NORMAL TIME---------------------------------------------------------------
        //Assumes time is normally distributed and does a kernel density estimation using
        //observation error variance as kernel width
        //Only works with RHF for now
        //return this.likelihoodFromNormalTime(time, ensembles, 100);
        //-------------------------DISCRETE TIME-------------------------------------------------------------
        //Creates a scaled kernel within each bin - approaches the correct value as bin quantity increases
        //Only works with RHF for now
        //return this.likelihoodFromDiscreteTime(time, ensembles);
        //-------------------------NORMAL LIKELIHOOD---------------------------------------------------------
        //Using the Monte Carlo method, generate some number of pseudoObservations
        //Then fit them all to a normal
        //Only works with EAKF or other Gaussian-assuming filters for now
        return this.normalLikelihood(time, ensembles, 100); 
        //------------------------KNOWN ERROR NORMAL---------------------------------------------------------
        //Applies the Normal Likelihood method, but does not attempt to infer likelihood, instead using known values
        //Only works with EAKF for now; requires minor changes to code for now
        //return this.knownErrorNormalLikelihood(time, ensembles, 100);
    }

    /** 
     * Fast alternative to normalLikelihood if error is known
     * Takes in observation time and a timeseries of the ensembles up until that time
     */
    Likelihood!dim knownErrorNormalLikelihood(double time, Timeseries!(Ensemble!dim) ensembles, uint kernels) {
        //Get the value of the observation at the given time
        Vector!(double, dim) obs = this.observations.value!dim(time);
        //These lists will store the values of the pseudo measurements independently
        //So that we can assimilate all variables separately
        double[][dim] pseudoMeasurements;
        //Do this a specified number of times: create a pseudo observation
        foreach(i; 0..kernels) {
            //Find a new time for the measurement
            auto newTime = NormalVariable!double(this.timeOffset, this.timeError);
            //Find a trajectory through the observation, then find the value on that trajectory
            //at the observation time (the pseudo-truth)
            Vector!(double, dim) base = this.integrator.integrateTo(obs, newTime(*this.gen), 1);
            NormalVariable!double normal;
            static foreach (j; 0..dim) {
                normal = NormalVariable!double(base[j], this.stateError[j]);
                pseudoMeasurements[j] ~= normal(*this.gen);
            }
        }
        //Find the mean of the pseudo-observations
        Vector!(double, dim) meanVector = new Vector!(double, dim)(pseudoMeasurements[].map!(a => a.mean).array);
        //Find the standard deviation of the pseudo-observations
        Vector!(double, dim) deviationVector = new Vector!(double, dim)(pseudoMeasurements[].map!(a => a.standardDeviation!1).array);
        //Return as a Gaussian likelihood with the above mean and standard deviation
        return new Likelihood!dim(meanVector, deviationVector);
    }

    /**
     * Returns a normal likelihood from normal-assumed time
     * TODO: Move into its own LikelihoodGetter
     */
    Likelihood!dim normalLikelihood(double time, Timeseries!(Ensemble!dim) ensembles, uint kernels) {
        //Get the most recent ensemble
        Ensemble!dim ensemble = new Ensemble!dim(ensembles.members[$ - 1].members);
        //If the ensemble is not yet past the maximum offset, integrate it through the interval:
        if(ensembles.times[$ - 1] + maximumOffset > ensembles.times[$ - 1]) {
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
        }
        //Make sure the observation time fits an observation
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        //Get the ensemble at the observation time
        ensemble = ensembles.value(time, this.integrator);
        //Infer the time likelihood
        this.getTimeLikelihood(time, ensembles);
        //Get the observation at the given time
        Vector!(double, dim) obs = this.observations.value!dim(time);
        //These lists will store the values of the pseudo measurements independently
        //So that we can assimilate all variables separately
        double[][dim] pseudoMeasurements;
        //Do this a specified number of times: create a pseudo observation
        foreach(i; 0..kernels) {
            //Find a new time for the measurement
            auto newTime = NormalVariable!double(this.timeOffset, this.timeDeviation);
            //Find a trajectory through the observation, then find the value on that trajectory
            //at the observation time (the pseudo-truth)
            Vector!(double, dim) base = this.integrator.integrateTo(obs, newTime(*this.gen), 1);
            NormalVariable!double normal;
            static foreach (j; 0..dim) {
                normal = NormalVariable!double(base[j], this.stateError[j]);
                pseudoMeasurements[j] ~= normal(*this.gen);
            }
        }
        //Find the mean of the pseudo-observations
        Vector!(double, dim) meanVector = new Vector!(double, dim)(pseudoMeasurements[].map!(a => a.mean).array);
        //Find the standard deviation of the pseudo-observations
        Vector!(double, dim) deviationVector = new Vector!(double, dim)(pseudoMeasurements[].map!(a => a.standardDeviation!1).array);
        //Return as a Gaussian likelihood with the above mean and standard deviation
        return new Likelihood!dim(meanVector, deviationVector);
    }

    /**
     * Returns likelihood packaged with discretely defined experimentally determined pdf for a given time
     * Fits time likelihood to a Gaussian, then uses a Monte Carlo-style kernel density approximation with kernel widths 
     * given by manufacturer
     * Assumes time likelihood is Gaussian (data suggest it will be with current method; it is advised to use likelihoodFromDiscrete time otherwise)
     */
    Likelihood!dim likelihoodFromNormalTime(double time, Timeseries!(Ensemble!dim) ensembles, uint kernels) {
        Ensemble!dim ensemble = new Ensemble!dim(ensembles.members[$ - 1].members);
        //If the ensemble is not yet past the maximum offset, integrate it through the interval:
        if(ensembles.times[$ - 1] + maximumOffset > ensembles.times[$ - 1]) {
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
        }
        //Ensure that the observation time has an associated observation
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        //Get the ensemble at the observation time
        ensemble = ensembles.value(time, this.integrator);
        //Infer time likelihood
        this.getTimeLikelihood(time, ensembles);
        //We define a discrete likelihood separately for the three variables 
        //in order to assimilate them separately
        double[][dim] likelihoods;
        //We give them each a miniscule value so that we don't end up with a uniformly zero likelihood
        //still, if this is happening, there is a problem...
        static foreach (i; 0..dim) {
            foreach(j; 0..ensemble.size) {
                likelihoods[i][j] = 0.000000001;
            }
        }
        //Get the observation at the give time
        Vector!(double, dim) obs = this.observations.value!dim(time);
        double[] pseudoTimes;
        //Generates pseudotimes; number of kernels can be changed within the code, but should probably be optimized
        //Generate some pseudo-observation times; these should be distributed normally
        foreach(i; 0..kernels) {
            auto newTimeVar = NormalVariable!double(this.timeOffset, this.timeDeviation);
            pseudoTimes ~= clamp(newTimeVar(*this.gen), this.minimumOffset, this.maximumOffset);
        }
        //We can create placeholder vectors here with the time as the x component
        //This way as the list is iterated through the x component will be overwritten by a location
        //and no extra memory assignment is necessary
        Vector!(double, dim)[] bases = pseudoTimes.dup.map!(a => new Vector!(double, dim)(a)).array;
        //In parallel, iterate through the placeholders
        foreach(ref component; bases.parallel) {
            //Fill the placeholder vector with the position of the obs trajectory at the given time
            component = this.integrator.integrateTo(obs, component[0], 1);
            //In parallel, iterate through the probabilities
            foreach(index, ref pointProbability; likelihoods[0].to!(double[]).parallel) {
                //Add the value of the normal probability density function with mean at the location being considered and error given a priori
                //to the probability in each bin
                static foreach (i; 0..dim) {
                    likelihoods[i][index] += normalVal(ensemble.members[index][i], component[i], this.stateError[i]);
                }
                //We end up with the sum of the normal pdfs centered at each center point in each bin
            }
        }
        //Normalize the likelihoods
        double[] sums = likelihoods[].map!(a => a.to!(double[]).sum).array;
        assert(sums.to!(double[]).all!(a => a != 0), "Experimental likelihood sum is 0");
        foreach(index, ref component; likelihoods[0].parallel ) {
            static foreach (i; 0..dim) {
                likelihoods[i][index] /= sums[i];
            }
        }
        //Return a discrete likelihood
        return new Likelihood!dim(likelihoods);
    }

    /**
     * Gets time likelihood, accounting for ensemble variance
     */
    void getTimeLikelihood(double time, Timeseries!(Ensemble!dim) ensembles) {
        //Ensure that the observation time is close enough to one that we know
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        //Get observation at time
        Vector!(double, dim) obs = this.observations.value!dim(time);
        //Split the interval into bins; TODO: Perhaps make this a property
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        //Get the times in the centers of each bin
        double[] binMiddles = iota(time + this.minimumOffset + binWidth / 2, time + this.maximumOffset, binWidth).array;
        //Get the system state at the middle of each bin
        //Ensembles passed into function should already have been advanced to the end of the interval
        Vector!(double, dim)[] binValues = binMiddles.map!(a => ensembles.meanSeries!dim.value!dim(a, this.integrator)).array;
        //Get the values of the ensemble points at the bin middles;
        //Might be faster to do this, then compute means from ensembles
        Ensemble!dim[] ensembleValues = binMiddles.map!(a => ensembles.value!dim(a, this.integrator)).array;
        //Initialize an array to store time likelihood - uniformly zero for now
        double[] binQuantities; 
        //Can maybe do this faster or more elegantly
        foreach(i; 0..this.bins) {
            binQuantities ~= 0;
        }
        //The observation error covariance is 0 outside the diagonal because we assume obs error is dimensionally independent
        Matrix!(double, dim, dim) obsErrorCovariance = new Matrix!(double, dim, dim)(0);
        static foreach (i; 0..dim) {
            obsErrorCovariance[i][i] = this.stateError[i] * this.stateError[i];
        }
        foreach(i; 0..binMiddles.length) {
            Vector!(double, dim) distance = binValues[i] - obs;
            //Be careful; options 1 and 2 are biased because they do not account for ensemble variance
            //However, they are computationally faster since option 3 involves matrix inversion and multiplication

            //Option 1: Discrete binary
            //If within 3 standard deviations, use it
            //binQuantities += (this.stateError * 3 - distance).any!(a => a < 0)? 0 : 1

            //Option 2: Smooth
            //Assume error is dimensionally independent, then find probability by multiplying the probabilities
            //binQuantities[i] = trivariateNormalVal(obs, binValues[i], obsErrorCovariance);

            //Option 3: Smooth multivariate
            //Assume error is dimensionally independent, then find probability
            //By finding the value of a multivariate (in this case trivariate) normal pdf at the value of distance
            //The covariance matrix of the ensemble has the variance across the diagonals and the covariance everywhere else
            Matrix!(double, dim, dim) ensembleCovariance = covariance!(dim, 1)(ensembleValues[i].valueLists);
            //The value of the probability is given by the trivariate normal PDF at the observation with mean ensemble mean
            //and covariance the sum of the ensemble covariance and the obs error covariance
            //This is the probability of taking the observation if the truth is taken from multivariate normal distribution
            //that is represented by the ensemble
            binQuantities[i] = multivariateNormalVal!dim(obs, binValues[i], obsErrorCovariance + ensembleCovariance);
        }
        //In parallel, update the time likelihood with the new inferred likelihood
        if(this.multiply) {
            //If the bin quantities are uniformly zero skip this
            if (binQuantities.sum == 0) return;
            //Fit Gaussian to inferred distribution
            double mean = PDFMean(binMiddles, binQuantities) - time;
            double variance = PDFVariance(binMiddles, binQuantities);
            double priorMean = this.timeOffset;
            double priorVariance = this.timeError;
            this.timeError = priorVariance == 0 || variance == 0 ? 0 :
                sqrt(1 / (1 / priorVariance + 1 / variance));
            this.timeOffset = priorVariance == 0 ? priorMean : variance == 0 ? mean : 
                this.timeError * (priorMean / priorVariance + mean / variance);
        } else {
            //If the bin quantities are uniformly zero skip this
            if (binQuantities.sum == 0) return;
            //Fit Gaussian to inferred distribution
            double mean = PDFMean(binMiddles, binQuantities) - time;
            double variance = PDFVariance(binMiddles, binQuantities);
            //writeln(binMiddles.map!(a => a - time).array);
            //writeln(binQuantities);
            //writeln("Mean ", mean);
            writeln("Variance ", variance);
            if (this.correctSlope) {
                variance = max(0, variance - 2 * this.stateError[0].pow(2) / (this.integrator.slope(ensembles.meanSeries!dim.value(time, this.integrator)).magnitude.pow(2))); 
                writeln("Correcting factor ", 2 * this.stateError[0].pow(2) / (this.integrator.slope(ensembles.meanSeries!dim.value(time, this.integrator)).magnitude.pow(2)));
                writeln("Corrected Variance ", variance);
            }
            this.timeError = (this.timeDeviation + sqrt(variance)).pow(2)/4;
            this.timeOffset = (this.timeOffset * this.stepsTaken + mean) / (this.stepsTaken + 1);
            //writeln("New inferred dist: N( ", this.timeOffset, ", ", this.timeError, " )");
            //File("data/sandbox/test.csv", "a").writeln(this.timeOffset, ", ", this.timeError);
        }
        this.stepsTaken += 1;
    }

    /**
     * Gets baselines for each bin to find likelihood for a given observation 
     * Backwards integration may cause instability; be careful!
     */
    Vector!(double, dim)[] getObservationTrajectory(double time) {
        //Ensure that there is an observation at the given time
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in observation times");
        //Find the observation at that time
        Vector!(double, dim) obs = this.observations.value!dim(time);
        Vector!(double, dim)[] baselines;
        const double binWidth = (this.maximumOffset - this.minimumOffset) / this.bins;
        //Find the value of the trajectory at the minimum offset time
        Vector!(double, dim) binEdge = this.integrator.integrateTo(obs, this.minimumOffset, 10 * cast(uint) (abs(this.minimumOffset) / binWidth));
        //Find the value of the trajectory at the middle of the first bin
        Vector!(double, dim) state = this.integrator.integrate(binEdge, binWidth / 2);
        baselines ~= state;
        //Find the values of the trajectory at every other bin by integrating forward with step binWidth
        foreach(i; 0..this.bins - 1) {
            state = this.integrator.integrate(state, binWidth);
            baselines ~= state;
        }
        return baselines;
    }

}