/**
 * A likelihood getter is a class that takes in information about the system and a time, then eturns a likelihood
 * This one is a parent class for all non-standard likelihood getters
 * It also acts as its own likelihood getter, which returns a standard Gaussian likelihood
 * This is the most basic form of likelihood available, and the most popular one
 * Perhaps this should be merged with DiscreteGaussianLikelihood, although DiscreteGaussianLikelihood is not currently in use
 */
module assimilation.likelihood.LikelihoodGetter;

import std.algorithm; //Used to check if an observation exists at the current time
import std.math; //Used to check approximate equality among doubles
import assimilation.likelihood.Likelihood; //Used to output a standard data storage class
import data.Ensemble; //Used to standardize discrete likelihood output, though this class does not use them
import data.Timeseries; //Used to keep track of observations

/**
 * Gets the likelihood for assimilation given information about the system
 */
class LikelihoodGetter {

    Timeseries!Vector observations; /// A list of all past observations that the likelihood getter is able to get likelihood for
    Vector stateError; //A vector representing the standard deviation of the observation error, which should be given a priori by the instrument manufacturer

    /**
     * Constructs the likelihood getter with a priori information about observations
     */ 
    this(Timeseries!Vector observations, Vector stateError) {
        this.observations = observations;
        this.stateError = stateError; //We assume this is constant over time, but it might not be
    }

    /**
     * Returns likelihood packaged with mean and deviation for a given time
     */
    Likelihood opCall(double time) {
        //Ensure an observation exists at that time
        assert(this.observations.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)), "Time not in list");
        return new Likelihood(this.observations.value(time), this.stateError);
    }

    /**
     * Returns likelihood given ensembles
     */
    Likelihood opCall(double time, Timeseries!Ensemble ensembles) {
        return null;
    }

}
