/**
 * This gets the values of a Gaussian likelihood at the ensemble members
 * It is useful for when assimilating with a normal likelihood and the Rank Histogram Filter
 * because RHF takes in a value of likelihood for each ensemble member
 */
module assimilation.likelihood.DiscreteGaussianLikelihood;

import std.algorithm; //Used for map function
import std.array; //Used to convert map output to array
import assimilation.likelihood.Likelihood; //Used to handle outputs from this class's opCall
import assimilation.likelihood.LikelihoodGetter; //Used for the parent class of this class
import data.Ensemble; //Used to define likelihood only at ensemble members
import data.Timeseries; //Used to standardize inputs across all LikelihoodGetter children
import integrators.Integrator; //Used to get values of otherwise undefined timeseries
import math.Vector; //Used for data storage
import utility.Normal; //Used to find the value of the normal probability density function for discrete likelihood

/**
 * A likelihood getter that returns a discrete representation of a Gaussian likelihood defined at the ensemble points for a given time
 * Useful for e.g. RHF
 * Knows system dimensions at compiletime
 */
class DiscreteGaussianLikelihood(uint dim) : LikelihoodGetter {

    Timeseries!(Ensemble!dim) ensembles; ///A timeseries of all past ensembles, in order to get the ensemble at the observation time
    Integrator integrator; ///A differential equation solver for finding ensemble values when they are not yet defined

    /**
     * Constructs the likelihood getter to prepare it for getting likelihood
     */
    this(Timeseries!(Vector!(double, dim)) observations, Timeseries!(Ensemble!dim) ensembles, Vector!(double, dim) stateError, Integrator integrator) {
        super(observations, stateError);
        this.ensembles = ensembles;
        this.integrator = integrator;
    }

    /**
     * Returns likelihood packaged with discretely defined Gaussian likelihoods for a given time
     */
    override Likelihood opCall(double time) {
        assert(this.observations.times.canFind(time));
        Ensemble!dim ensemble = ensembles.value(time, this.integrator);
        //Find the value of the Gaussial likelihood pdf at each ensemble point
        double[][] likelihoods;
        static foreach (i; 0..dim) {
            likelihoods[i] = ensemble.valueLists[i].map!(a => normalVal(a, this.observations.timeAssociate[time][i], this.stateError[i])).array;
        }
        return new Likelihood(likelihoods);
    }

}