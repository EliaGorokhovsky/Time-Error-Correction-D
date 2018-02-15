module experiment.Experiment;

import std.algorithm;
import std.datetime;
import std.math;
import std.range;
import std.typecons;
import std.file;
import assimilation.Assimilator;
import assimilation.EAKF;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import data.Ensemble;
import data.Timeseries;
import data.Vector;
import experiment.Analytics;
import experiment.error.ErrorGenerator;
import integrators.Integrator;
import systems.System;

/**
 * A class responsible for handling the experiment itself
 */
class Experiment {

    Integrator integrator;
    ErrorGenerator errorGen;
    Assimilator assimilator;
    LikelihoodGetter likelihoodGetter;

    Timeseries!Vector truth;
    Timeseries!Vector observations;
    Timeseries!Ensemble ensembleSeries;

    this(Integrator integrator, Assimilator assimilator) {
        this.integrator = integrator;
        this.assimilator = assimilator;
    }

    /** 
     * Sets the likelihood getter for the experiment
     * Must be done after getTruth and getObservations
     */
    void setLikelihood(LikelihoodGetter likelihoodGetter) {
        this.likelihoodGetter = likelihoodGetter;
    }

    /** 
     * Sets the error for the experiment
     */
    void setError(ErrorGenerator errorGen) {
        this.errorGen = errorGen;
    }

    /**
     * Integrates a start point over the specified time interval
     * Ends at the last timestep before endTime
     */
    Timeseries!Vector getTruth(Vector point, double startTime, double endTime, double dt) {
        Timeseries!Vector truth = new Timeseries!Vector([point], [0]);
        foreach(i; iota(startTime, endTime, dt)) {
            point = this.integrator(point, dt);
            truth.add(i + dt, point);
        }
        this.truth = truth;
        return this.truth;
    }

    /**
     * Gets a set of observations 
     * Ends at the last interval before endTime
     */
    Timeseries!Vector getObservations(double startTime, double endTime, double interval) {
        Timeseries!Vector observations = new Timeseries!Vector();
        foreach(i; iota(startTime, endTime, interval)) {
            observations.add(i, this.errorGen(i));
        }
        this.observations = observations;
        return this.observations;
    }

    /**
     * Runs an ensemble over an interval
     * Spin-up time: no assimilation
     */
    Timeseries!Ensemble getEnsembleTimeseries(bool experiment)(double startTime, double endTime, double dt, double spinup, Ensemble ensemble) {
        Timeseries!Ensemble ensembleSeries = new Timeseries!Ensemble();
        ensembleSeries.add(0, ensemble);
        assert(ensembleSeries.members !is null, "Ensemble series is null");
        foreach(i; iota(startTime, endTime, dt)) {
            if(this.observations.times.any!(a => a.approxEqual(i, 1e-06, 1e-06)) && i >= spinup) {
                ensemble *= 1.5;
                Timeseries!Ensemble placeholder = new Timeseries!Ensemble(ensembleSeries.members, ensembleSeries.times);
                this.assimilator.setLikelihood(experiment? this.likelihoodGetter(i, placeholder) : this.likelihoodGetter(i));
                ensemble = this.assimilator(ensemble);
            }
            ensemble = this.integrator(ensemble, dt);
            ensembleSeries.add(i + dt, ensemble);
        }
        this.ensembleSeries = ensembleSeries;
        return this.ensembleSeries;
    }

}

unittest {

    import std.stdio;
    import integrators.RK4;
    import experiment.Analytics;
    import experiment.error.GaussianError;
    import assimilation.EAKF;

    writeln("\nUNITTEST: Experiment");
    class Test : System {
        override Vector opCall(Vector state) { return Vector(1, 1, 1); }
    }
    RK4 rk4 = new RK4(new Test());
    Experiment process = new Experiment(rk4, new EAKF);
    process.getTruth(Vector(0, 0, 0), 0, 10, 1);
    writeln("Integrating <1, 1, 1> from 0 to 10 returns ", process.truth.members);
    process.setError(new GaussianError(Vector(0.1, 0.1, 0.1), process.truth, rk4));
    process.getObservations(0, 10, 3);
    writeln("Observing every 3 seconds with std (0.1, 0.1, 0.1) returns ", process.observations.members);
    process.setLikelihood(new LikelihoodGetter(process.observations, Vector(0.1, 0.1, 0.1)));
    process.getEnsembleTimeseries!false(0, 10, 1, 4, new Ensemble(Vector(0, 0, 0), 3, Vector(0.1, 0.1, 0.1)));
    writeln("RMSE is ", RMSE(process.ensembleSeries, process.truth));

}