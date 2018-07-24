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
import utility.Normal; //Used to find the value of the normal probability density function for discrete likelihood

/**
 * A likelihood getter that returns a discrete representation of a Gaussian likelihood defined at the ensemble points for a given time
 * Useful for e.g. RHF
 */
class DiscreteGaussianLikelihood : LikelihoodGetter {

    Timeseries!Ensemble ensembles; ///A timeseries of all past ensembles, in order to get the ensemble at the observation time
    Integrator integrator; ///A differential equation solver for finding ensemble values when they are not yet defined

    /**
     * Constructs the likelihood getter to prepare it for getting likelihood
     */
    this(Timeseries!Vector observations, Timeseries!Ensemble ensembles, Vector stateError, Integrator integrator) {
        super(observations, stateError);
        this.ensembles = ensembles;
        this.integrator = integrator;
    }

    /**
     * Returns likelihood packaged with discretely defined Gaussian likelihoods for a given time
     */
    override Likelihood opCall(double time) {
        assert(this.observations.times.canFind(time));
        Ensemble ensemble = ensembles.value(time, this.integrator);
        //Find the value of the Gaussial likelihood pdf at each ensemble point
        return new Likelihood(
            ensemble.xValues.map!(a => normalVal(a, this.observations.timeAssociate[time].x, this.stateError.x)).array,
            ensemble.yValues.map!(a => normalVal(a, this.observations.timeAssociate[time].y, this.stateError.y)).array,
            ensemble.zValues.map!(a => normalVal(a, this.observations.timeAssociate[time].z, this.stateError.z)).array
        );
    }

}