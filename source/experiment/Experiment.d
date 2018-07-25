module experiment.Experiment;

import std.algorithm;
import std.datetime;
import std.math;
import std.range;
import std.stdio;
import std.typecons;
import std.file;
import assimilation.Assimilator;
import assimilation.EAKF;
import assimilation.likelihood.DiscreteExperimentalLikelihood;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import data.Ensemble;
import data.Timeseries;
import experiment.Analytics;
import experiment.error.ErrorGenerator;
import integrators.Integrator;
import math.Vector;
import systems.System;

/**
 * A class responsible for handling the experiment itself
 */
class Experiment(uint dim) {

    Integrator integrator;
    ErrorGenerator!dim errorGen;
    Assimilator!dim assimilator;
    EAKF!dim standardEAKF = new EAKF();
    LikelihoodGetter!dim standardLikelihood;
    LikelihoodGetter!dim likelihoodGetter;

    Timeseries!(Vector!(double, dim)) truth;
    Timeseries!(Vector!(double, dim)) observations;
    Timeseries!(Ensemble!dim) ensembleSeries;

    this(Integrator integrator, Assimilator!dim assimilator) {
        this.integrator = integrator;
        this.assimilator = assimilator;
    }

    /** 
     * Sets the likelihood getter for the experiment
     * Must be done after getTruth and getObservations
     */
    void setLikelihood(LikelihoodGetter!dim likelihoodGetter) {
        this.likelihoodGetter = likelihoodGetter;
    }

    /** 
     * Sets the error for the experiment
     */
    void setError(ErrorGenerator!dim errorGen) {
        this.errorGen = errorGen;
    }

    /**
     * Integrates a start point over the specified time interval
     * Ends at the last timestep before endTime
     */
    Timeseries!(Vector!(double, dim)) getTruth(Vector!(double, dim) point, double startTime, double endTime, double dt) {
        Timeseries!(Vector!(double, dim)) truth = new Timeseries!(Vector!(double, dim))([point], [0]);
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
    Timeseries!(Vector!(double, dim)) getObservations(double startTime, double endTime, double interval) {
        Timeseries!(Vector!(double, dim)) observations = new Timeseries!(Vector!(double, dim))();
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
    Timeseries!(Ensemble!dim) getEnsembleTimeseries(bool experiment)(double startTime, double endTime, double dt, double spinup, double priming, Ensemble!dim ensemble) {
        Timeseries!(Ensemble!dim) ensembleSeries = new Timeseries!(Ensemble!dim)();
        ensembleSeries.add(0, ensemble);
        assert(ensembleSeries.members !is null, "Ensemble series is null");
        foreach(i; iota(startTime, endTime, dt)) {
            if(this.observations.times.any!(a => a.approxEqual(i, 1e-06, 1e-06)) && i >= spinup) {
                //ensemble *= 1.5;
                Timeseries!(Ensemble!dim) placeholder = new Timeseries!(Ensemble!dim)(ensembleSeries.members, ensembleSeries.times);
                this.assimilator.setLikelihood(experiment? this.likelihoodGetter(i, placeholder) : this.likelihoodGetter(i));
                if(experiment && i < priming) {
                    this.standardEAKF.setLikelihood(this.standardLikelihood(i));
                    ensemble = this.standardEAKF(ensemble);
                }
                else {
                    ensemble = this.assimilator(ensemble);
                }
            }
            /*if(i % 5 == 0) { 
                writeln("Time ", i, " for ", experiment? "treatment" : "control");
            }*/
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
    import mir.random;

    writeln("\nUNITTEST: Experiment");
    class Test : System {
        override Vector!(double, 3) opCall(Vector!(double, 3) state) { return new Vector!(double, 3)(1); }
    }
    RK4 rk4 = new RK4(new Test());
    Experiment!3 process = new Experiment!3(rk4, new EAKF!3);
    process.getTruth(new Vector!(double, 3)(0), 0, 10, 1);
    writeln("Integrating <1, 1, 1> from 0 to 10 returns ", process.truth.members);
    Random gen = Random(unpredictableSeed);
    process.setError(new GaussianError!3(new Vector!(double, 3)(0.1), process.truth, rk4, &gen));
    process.getObservations(0, 10, 3);
    writeln("Observing every 3 seconds with std (0.1, 0.1, 0.1) returns ", process.observations.members);
    process.setLikelihood(new LikelihoodGetter(process.observations, new Vector!(double, 3)(0.1)));
    process.getEnsembleTimeseries!false(0, 10, 1, 4, 0, new Ensemble(new Vector!(double, 3)(0), 3, new Vector!(double, 3)(0.1), &gen));
    writeln("RMSE is ", RMSE(process.ensembleSeries, process.truth));

}