module assimilation.LikelihoodGetter;

import std.algorithm;
import std.typecons;
import assimilation.Likelihood;
import data.Ensemble;
import data.Timeseries;
import data.Vector;

/**
 * Gets the likelihood for assimilation given information about the system
 */
class LikelihoodGetter {

    Timeseries!Vector observations;
    Vector stateError;
    double timeError;

    this(Timeseries!Vector observations, Vector stateError, double timeError) {
        this.observations = observations;
        this.stateError = stateError;
        this.timeError = timeError;
    }

    /**
     * Returns observation and likelihood for a given time
     * Output for EAKF 
     */
    Likelihood opCall(double time) {
        assert(this.observations.times.canFind(time));
        return new Likelihood(this.observations.timeAssociate[time], this.stateError);
    }

}