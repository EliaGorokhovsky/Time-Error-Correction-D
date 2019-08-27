module experiment.Experiment;

import std.algorithm;
import std.datetime;
import std.math;
import std.range;
import std.typecons;
import std.file;
import std.stdio;
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
import math.Matrix;
import math.Vector;
import systems.System;

/**
 * A class responsible for handling the experiment itself
 */
class Experiment(uint dim) {

    Integrator!dim integrator;
    ErrorGenerator!dim errorGen;
    Assimilator!dim assimilator;
    EAKF!dim standardEAKF = new EAKF!dim();
    LikelihoodGetter!dim standardLikelihood;
    LikelihoodGetter!dim likelihoodGetter;

    Timeseries!(Vector!(double, dim)) truth;
    Timeseries!(Vector!(double, dim)) observations;
    Timeseries!double observationTimes;
    Timeseries!double inferredObservationTimes;
    Timeseries!(Ensemble!dim) ensembleSeries;

    /**
     * Returns the RMSE of the ensemble
     */
    @property double ensembleRMSE() {
        return RMSE!dim(this.ensembleSeries, this.truth);
    }

    this(Integrator!dim integrator, Assimilator!dim assimilator) {
        this.integrator = integrator;
        this.assimilator = assimilator;
    }

    /**
     * Returns the RMSE of the observation time reports
     */
    double getTimeErrorRMSE(ulong start) {
        return sqrt(
            iota(start, this.observationTimes.length, 1)
                .fold!((sum, i) => sum + pow(this.observationTimes.times[i] - this.observationTimes.members[i], 2))(0.0)
            / (this.observationTimes.length - start)
        );
    }

    /**
     * Returns the RMSE of the calculated time reports
     */
    double getInferredTimeErrorRMSE(ulong start) {
        return sqrt(
            iota(start, this.inferredObservationTimes.length, 1)
                .fold!((sum, i) => sum + pow(this.inferredObservationTimes.times[i] - this.observationTimes.valueAtTime(this.inferredObservationTimes.times[i]) + this.inferredObservationTimes.members[i], 2))(0.0)
            / (this.inferredObservationTimes.length - start)
        );
    }

    /**
     * Returns how frequently the algorithm correctly guesses whether the observation time is over- or under- estimated
     */
    double getDirectionGuessRate(ulong start) {
        return cast(double)(iota(start, this.inferredObservationTimes.length, 1)
                .count!(i => 
                (this.observationTimes.members[i + 1] < this.observationTimes.times[i + 1] && this.inferredObservationTimes.members[i] < 0) 
                || 
                (this.observationTimes.members[i + 1] > this.observationTimes.times[i + 1] && this.inferredObservationTimes.members[i] > 0)
            )()) / (this.inferredObservationTimes.length - start);
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
     * Sets the true state to the input
     */
    void setTruth(Timeseries!(Vector!(double, dim)) truth) {
        this.truth = truth;
    }

    /**
     * Gets a set of observations 
     * Ends at the last interval before endTime
     */
    Timeseries!(Vector!(double, dim)) getObservations(double startTime, double endTime, double interval) {
        Timeseries!(Vector!(double, dim)) observations = new Timeseries!(Vector!(double, dim))();
        Timeseries!double observationTimes = new Timeseries!double();
        foreach(i; iota(startTime, endTime, interval)) {
            Tuple!(double, Vector!(double, dim)) obs = this.errorGen.generate(i);
            observationTimes.add(i, obs[0]);
            observations.add(i, obs[1]);
        }
        this.observationTimes = observationTimes;
        this.observations = observations;
        return this.observations;
    }

    /**
     * Sets the observations for this run to the inputted ones
     */
    void setObservations(Timeseries!(Vector!(double, dim)) observations, Timeseries!double observationTimes) {
        this.observations = observations;
        this.observationTimes = observationTimes;
    }

    /**
     * Runs an ensemble over an interval
     * Spin-up time: no assimilation
     */
    Timeseries!(Ensemble!dim) getEnsembleTimeseries(bool experiment)(string testfilename, double startTime, double endTime, double dt, double spinup, double priming, Ensemble!dim ensemble) {  
        //File(testfilename, "a").writeln("Time, Observed time error, Predicted time error, Time error uncertainty, Truth,,, Observation,,, Ensemble Mean,,, Ensemble Variance,,, Likelihood Standard Deviation,,, Slope,,, Predicted time err var, Errsum, DiffErrsum, Inverse square time error, Time offset calculation denom");
        Timeseries!(Ensemble!dim) ensembleSeries = new Timeseries!(Ensemble!dim)();
        Timeseries!double inferredObservationTimes = new Timeseries!double();
        ensembleSeries.add(0, ensemble);
        assert(ensembleSeries.members !is null, "Ensemble series is null");
        foreach(i; iota(startTime, endTime, dt)) {
            if(this.observations.times.any!(a => a.approxEqual(i, 1e-06, 1e-06)) && i >= spinup) {
                //Prior inflation goes here.
                //ensemble *= 1.5;
                Timeseries!(Ensemble!dim) placeholder = new Timeseries!(Ensemble!dim)(ensembleSeries.members, ensembleSeries.times);
                //Sets the likelihood used by assimilation. Likelihood inflation goes here.
                Likelihood!dim likelihood = experiment? this.likelihoodGetter(i, placeholder) : new Likelihood!dim(this.likelihoodGetter(i).gaussianMean, this.likelihoodGetter(i).gaussianDeviation);
                //Likelihood!dim inflatedLikelihood = new Likelihood!dim(likelihood.gaussianMean, likelihood.gaussianDeviation * sqrt(2.0));
                this.assimilator.setLikelihood(likelihood);
                if(experiment && i < priming) {
                    this.standardEAKF.setLikelihood(this.standardLikelihood(i));
                    ensemble = this.standardEAKF(ensemble);
                }
                else {
                    ensemble = this.assimilator(ensemble);
                }
                //Record observation time and expected observation time.
                if (DiscreteExperimentalLikelihood!3 lik = cast(DiscreteExperimentalLikelihood!3) (this.likelihoodGetter)) {
                    Vector!(double, dim) obs = this.observations.valueAtTime(i);
                    Vector!(double, dim) slope = /*new Ensemble!dim(ensemble.members.map!(a => this.integrator.slope(a)).array).eMean;*/
                                                this.integrator.slope(ensemble.eMean);
                    Vector!(double, dim) errors = this.likelihoodGetter.stateError;
                    Vector!(double, dim) diff = obs - ensemble.eMean;
                    double errSum = iota(0, dim, 1).fold!((sum, i) => sum + slope[i] * slope[i] * pow(errors[i], -2))(0.0);
                    double diffErrSum = iota(0, dim, 1).fold!((sum, i) => sum + diff[i] * slope[i] * pow(errors[i], -2))(0.0);
                    //For unknown error
                    double timeError = lik.timeVariance/* - (dim + 1) / errSum*/;
                    //For known error
                    //double timeError = lik.timeError * lik.timeError;
                    double inverseSquareTimeError = 1 / timeError;
                    double denom = errSum + inverseSquareTimeError;
                    double mean = diffErrSum / denom;
                    double var = 1 / denom;

                    //if (abs(mean) <= 5000000 * sqrt(timeError)) {
                    inferredObservationTimes.add(i, mean);
                    //} else {*/
                        //inferredObservationTimes.add(i, this.observationTimes.valueAtTime(i) - i);
                    //}
                    /*File(testfilename, "a").writeln(
                        i, ", ", 
                        this.observationTimes.valueAtTime(i) - i, ", ", 
                        mean, ", ", 
                        var, ", ", 
                        this.truth.valueAtTime(i), ", ", 
                        obs, ", ", 
                        ensemble.eMean, ", ", 
                        ensemble.eVariance, ", ", 
                        likelihood.gaussianDeviation, ", ", 
                        slope, ", ", 
                        lik.timeVariance, ", ", 
                        errSum, ", ", 
                        diffErrSum, ", ", 
                        inverseSquareTimeError, ",", 
                        denom
                    );*/
                } 
            }
            /*if(i % 5 == 0) { 
                writeln("Time ", i, " for ", experiment? "treatment" : "control");
            }*/
            ensemble = this.integrator(ensemble, dt);
            ensembleSeries.add(i + dt, ensemble);
        }
        this.ensembleSeries = ensembleSeries;
        this.inferredObservationTimes = inferredObservationTimes;
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
    class Test : System!3 {
        override Vector!(double, 3) opCall(Vector!(double, 3) state) { return new Vector!(double, 3)(1); }
    }
    RK4!3 rk4 = new RK4!3(new Test());
    Experiment!3 process = new Experiment!3(rk4, new EAKF!3);
    process.getTruth(new Vector!(double, 3)(0), 0, 10, 1);
    writeln("Integrating <1, 1, 1> from 0 to 10 returns ", process.truth.members);
    Random gen = Random(1);
    process.setError(new GaussianError!3(new Vector!(double, 3)(0.1), process.truth, rk4, &gen));
    process.getObservations(0, 10, 3);
    writeln("Observing every 3 seconds with std (0.1, 0.1, 0.1) returns ", process.observations.members);
    process.setLikelihood(new LikelihoodGetter!3(process.observations, new Vector!(double, 3)(0.1)));
    Ensemble!3 base = new Ensemble!3(new Vector!(double, 3)(0), 3, new Vector!(double, 3)(0.1), &gen);
    writeln("Starting experiment with ensemble: ", base.members);
    process.getEnsembleTimeseries!false(0, 10, 1, 4, 0, base);
    writeln("RMSE is ", RMSE(process.ensembleSeries, process.truth));

}