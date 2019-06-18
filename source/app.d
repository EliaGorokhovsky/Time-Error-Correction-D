/*
All code written by Elia Gorokhovsky
*/
import std.algorithm;
import std.array;
import std.conv;
import std.datetime.stopwatch;
import std.file;
import std.math;
import std.range;
import std.stdio;
import mir.random;
import mir.random.variable;
import assimilation.Assimilator;
import assimilation.EAKF;
import assimilation.RHF;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import assimilation.likelihood.DiscreteExperimentalLikelihood;
import data.Ensemble;
import experiment.Analytics;
import experiment.Experiment;
import experiment.Parameters;
import experiment.error.GaussianError;
import experiment.error.GaussianTimeError;
import experiment.error.UniformTimeError;
import integrators.Integrator;
import integrators.RK4;
import math.Vector;
import systems.Line;
import systems.System;
import systems.Lorenz63;

enum dimensions = 3; ///How many dimensions the experiment is being run in
enum verboseRun = false; ///Whether to mention in detail what's going on in the program run

void run(Parameters!dimensions params, double observationInterval, double timeError, Vector!(double, dimensions) error, Random gen, ulong seed) {
	
	Vector!(double, dimensions) actualError = error; ///The standard deviation of the Gaussian error in space
	Vector!(double, dimensions) expectedError = error; ///The a priori expected standard deviation for Gaussian space error
	if(params.config == RunConfigurations.COMPARE_RMSE) {
		writeln("Control:");
		Experiment!dimensions control = new Experiment!dimensions(params.integrator, params.controlAssimilator);
		static if (verboseRun) writeln("Successfully initialized experiment.");
		//Get truth
		control.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		static if (verboseRun) writeln("Successfully got true timeseries.");
		//Get observations
		control.setError(new GaussianTimeError!dimensions(timeError, actualError, control.truth, params.integrator, &gen));
		static if (verboseRun) writeln("Successfully set error generator.");
		control.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		static if (verboseRun) writeln("Successfully generated observations.");
		//Do assimilation
		control.setLikelihood(new LikelihoodGetter!dimensions(control.observations, expectedError));
		static if (verboseRun) writeln("Successfully set likelihood.");
		control.getEnsembleTimeseries!false(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double controlRMSE = RMSE!dimensions(control.ensembleSeries, control.truth);
		writeln("Control RMSE for time error ", timeError, " is ", controlRMSE);

		writeln("Experiment:");
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		static if (verboseRun) writeln("Successfully initialized experiment.");
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		static if (verboseRun) writeln("Successfully got true timeseries.");
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		static if (verboseRun) writeln("Successfully set error generator.");
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		static if (verboseRun) writeln("Successfully generated observations.");
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood!dimensions(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen, //0, timeError
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter!dimensions(treatment.observations, expectedError);
		static if (verboseRun) writeln("Successfully set likelihood.");
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ",", error[0], ", ", timeError, ", ", controlRMSE, ", ", treatmentRMSE);
		writeln("Treatment RMSE for time error ", timeError, " is ", treatmentRMSE);
	} else if(params.config == RunConfigurations.CONTROL_RMSE) {
		writeln("Control:");
		Experiment!dimensions control = new Experiment!dimensions(params.integrator, params.controlAssimilator);
		static if (verboseRun) writeln("Successfully initialized experiment.");
		//Get truth
		control.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		static if (verboseRun) writeln("Successfully got true timeseries.");
		//Get observations
		control.setError(new GaussianTimeError!dimensions(timeError, actualError, control.truth, params.integrator, &gen));
		static if (verboseRun) writeln("Successfully set error generator.");
		control.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		static if (verboseRun) writeln("Successfully generated observations.");
		//Do assimilation
		control.setLikelihood(new LikelihoodGetter!dimensions(control.observations, expectedError));
		static if (verboseRun) writeln("Successfully set likelihood.");
		control.getEnsembleTimeseries!false(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double controlRMSE = RMSE!dimensions(control.ensembleSeries, control.truth);
		writeln("Control RMSE for time error ", timeError, " is ", controlRMSE);
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ",", error[0],  ", ", timeError, ", ", controlRMSE);
	} else if(params.config == RunConfigurations.TREATMENT_RMSE) {
		writeln("Experiment:");
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		static if (verboseRun) writeln("Successfully initialized experiment.");
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		static if (verboseRun) writeln("Successfully got true timeseries.");
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		static if (verboseRun) writeln("Successfully set error generator.");
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		static if (verboseRun) writeln("Successfully generated observations.");
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood!dimensions(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter!dimensions(treatment.observations, expectedError);
		static if (verboseRun) writeln("Successfully set likelihood.");
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ",", error[0], ", ", timeError, ", " , treatmentRMSE, ", ");
		writeln("Treatment RMSE for time error ", timeError, " is ", treatmentRMSE);
	} else if(params.config == RunConfigurations.INFERRED_TIME_ERROR) {
		writeln("Experiment:");
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		static if (verboseRun) writeln("Successfully initialized experiment.");
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		static if (verboseRun) writeln("Successfully got true timeseries.");
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		static if (verboseRun) writeln("Successfully set error generator.");
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		static if (verboseRun) writeln("Successfully generated observations.");
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood!dimensions(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter!dimensions(treatment.observations, expectedError);
		static if (verboseRun) writeln("Successfully set likelihood.");
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 2, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		DiscreteExperimentalLikelihood!dimensions treatmentLikelihood = cast(DiscreteExperimentalLikelihood!dimensions) treatment.likelihoodGetter;
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ",", error[0], ", ", timeError, ", ", treatmentLikelihood.expectedTime, ",", treatmentLikelihood.timeDeviation, ",", treatment.getTimeErrorRMSE(100), ",", treatment.getInferredTimeErrorRMSE(100), ",", treatment.getDirectionGuessRate(100), ",,", treatmentLikelihood.timeLikelihood.to!string[1 .. $ - 1]);
		writeln("Inferred time error for time error ", timeError, " is ", treatmentLikelihood.timeDeviation, " with observation time RMSE ", treatment.getTimeErrorRMSE(100), " and RMSE of calculated time error ", treatment.getInferredTimeErrorRMSE(100), " and direction guess rate ", treatment.getDirectionGuessRate(100));
	} else if (params.config == RunConfigurations.AVERAGE_DISTANCE) {
		writeln("Experiment:");
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, params.obsEndTime - params.obsStartTime + 10);
		treatment.setLikelihood(new LikelihoodGetter!dimensions(treatment.observations, expectedError));
		treatment.getEnsembleTimeseries!false(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.ensembleEndTime, params.ensembleEndTime, new Ensemble!dimensions(params.ensembleGenesis, 5, params.ensembleDeviation, &gen)
		);
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		writeln("Treatment RMSE: ", treatmentRMSE);
	}

}

enum RunConfigurations: string {
	COMPARE_RMSE = "Compare RMSE",
	CONTROL_RMSE = "Control RMSE",
	TREATMENT_RMSE = "Treatment RMSE",
	INFERRED_TIME_ERROR = "Inferred Time Error",
	AVERAGE_DISTANCE = "Average Distance"
}

void main() {
	RunConfigurations config = RunConfigurations.INFERRED_TIME_ERROR;
	string filename = "data/datacollection/Inferred-Time-Error-13-Electric-Boogaloo.csv";
	string logfile = "data/ExperimentLog.txt";
	string tag = "Testing for outliers";
	StopWatch stopwatch = StopWatch(AutoStart.no);
	bool logThisExperiment = true; //Set this to false if you don't want to write the experiment to the logfile
	double[] observationIntervals = [0.1];
	double[] timeErrors = [0.014];
	/*foreach(i; 0..15) {
		timeErrors ~= i * 0.001;
	}*/
	double[] errors = [0.1];
	/*foreach(i; 1..500) {
		errors ~= i * 0.02;
	}*/
	//This will set up a number of random seeds
	//The first map statement will give different random seeds every program run
	//The second map statement will ensure that all program runs are the same
	//You can also set random seeds to those outputted by the program to replicate its results
	ulong[] seeds = /*iota(0, 10, 1)*/
					/*.map!(a => unpredictableSeed)*/
					/*.map!(a => cast(ulong)a)*/
					/*.array;*/
					[8370133083869323966uL]
					/*[82029386284230530uL, 13963964892208834654uL, 9135740542362501819uL, 7031332123152652807uL, 15022117507543414108uL, 6622178055635864012uL, 1522060868933556120uL, 16100634437019494802uL, 8370133083869323966uL, 8296940596219331729uL]*/;
	//Package the parameters into one object
	Parameters!dimensions params = Parameters!dimensions(
		new Vector!(double, dimensions)(1), //The initial point of the truth
		0, //The initial time with which to associate the initial point
		80, //The time at which to stop the experiment
		0.01, //The length of each step of the integrator
		new RK4!dimensions(new Lorenz63()), //The integrator used to return points from previous points, and its system
		0, //When to start observing
		80, //When to stop observing
		0, //When to create the ensemble
		80, //When to stop assimilating
		0.01, //The step for ensemble integration
		0.1, //The amount of time the ensemble is run before beginning to assimilate
		new Vector!(double, dimensions)(1), //The mean of the initial ensemble distribution
		new Vector!(double, dimensions)(0.1), //The standard deviation of the initial ensemble distribution
		20, //The size of the ensemble
		new EAKF!dimensions(), //The assimilation method for the control
		new EAKF!dimensions(), //The assimilation method for the treatment 
		-0.06, //The first time that is a valid time for observation relative to reported time
		0.06, //The last time that is a valid time for observation relative to reported time
		11, //The amount of different time intervals tested in experimental likelihood algorithm
		observationIntervals, //The intervals between observations that will be tested
		timeErrors, //The time error standard deviations that will be tested
		errors, //The observation error standard deviations that will be tested
		seeds, //The seeds for each of the trials
		filename, //The file to write the data to
		config, //The run configuration
		tag //User-specified experiment information
	);
	writeln(params);
	if(logThisExperiment) File(logfile, "a").writeln(params);
	stopwatch.start();
	if(config == RunConfigurations.COMPARE_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Control RMSE, Treatment RMSE");
	} else if(config == RunConfigurations.CONTROL_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Control RMSE");
	} else if(config == RunConfigurations.TREATMENT_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Treatment RMSE");
	} else if(config == RunConfigurations.INFERRED_TIME_ERROR) {
		string binMiddles = iota(0, params.bins, 1).map!(a => params.minimumOffset + (params.maximumOffset - params.minimumOffset) / (params.bins * 2) + a * (params.maximumOffset - params.minimumOffset) / (params.bins)).array.to!string[1 .. $ - 1];
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Inferred Time Offset, Inferred Time Error, Observation Time RMSE, Adjusted Observation Time RMSE, Direction Guess Frequency, Inferred Time Offset Distribution: , ", binMiddles);
	} else if(config == RunConfigurations.AVERAGE_DISTANCE){}
	uint counter = 0;
	foreach(observationInterval; observationIntervals) {
		foreach(timeError; timeErrors) {
			foreach(error; errors) {
				foreach(ref seed; seeds) {
					run(params, observationInterval, timeError, new Vector!(double, dimensions)(error), Random(seed), seed);
					counter++;
					writeln("Experiment progress: ", counter, "/", observationIntervals.length * timeErrors.length * errors.length * seeds.length, " runs = ", 100 * counter / (observationIntervals.length * timeErrors.length * errors.length * seeds.length), "%");
				}
			}
		}
	}
	stopwatch.stop();
	writeln("Finished ", 
			observationIntervals.length * timeErrors.length * errors.length * seeds.length, 
			" trials in ", stopwatch.peek.total!"msecs", " ms.");
	if(logThisExperiment) File(logfile, "a").writeln("Complete in " ~ stopwatch.peek.total!"msecs".to!string ~ " ms!\n");
}