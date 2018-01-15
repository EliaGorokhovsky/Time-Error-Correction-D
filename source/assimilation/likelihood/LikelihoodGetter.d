module assimilation.likelihood.LikelihoodGetter;

import std.algorithm;
import assimilation.likelihood.Likelihood;
import data.Ensemble;
import data.Timeseries;
import data.Vector;

/**
 * Gets the likelihood for assimilation given information about the system
 */
class LikelihoodGetter {

    Timeseries!Vector observations;
    Vector stateError;

    this(Timeseries!Vector observations, Vector stateError) {
        this.observations = observations;
        this.stateError = stateError;
    }

    /**
     * Returns likelihood packaged with mean and deviation for a given time
     */
    Likelihood opCall(double time) {
        assert(this.observations.times.canFind(time));
        return new Likelihood(this.observations.timeAssociate[time], this.stateError);
    }

}
