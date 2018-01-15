module assimilation.likelihood.DiscreteGaussianLikelihood;

import std.algorithm;
import std.array;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import data.Ensemble;
import data.Timeseries;
import data.Vector;
import integrators.Integrator;
import utility.Normal;

/**
 * A likelihood getter that returns a discrete representation of a Gaussian likelihood defined at the ensemble points for a given time
 * Useful for e.g. RHF
 */
class DiscreteGaussianLikelihood : LikelihoodGetter {

    Timeseries!Ensemble ensembles;
    Integrator integrator;

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
        return new Likelihood(
            ensemble.xValues.map!(a => normalVal(a, this.observations.timeAssociate[time].x, this.stateError.x)).array,
            ensemble.yValues.map!(a => normalVal(a, this.observations.timeAssociate[time].y, this.stateError.y)).array,
            ensemble.zValues.map!(a => normalVal(a, this.observations.timeAssociate[time].z, this.stateError.z)).array
        );
    }

}