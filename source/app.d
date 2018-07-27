import std.algorithm;
import std.array;
import std.conv;
import std.datetime.stopwatch;
import std.file;
import std.math;
import std.parallelism;
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

enum dimensions = 1; ///How many dimensions the experiment is being run in

string run(Parameters!dimensions params, double observationInterval, double timeError, Vector!(double, dimensions) error, Random gen, ulong seed) {
	
	Vector!(double, dimensions) actualError = error; ///The standard deviation of the Gaussian error in space
	Vector!(double, dimensions) expectedError = error; ///The a priori expected standard deviation for Gaussian space error
	if(params.config == RunConfigurations.COMPARE_RMSE) {
		Experiment!dimensions control = new Experiment!dimensions(params.integrator, params.controlAssimilator);
		//Get truth
		control.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		//Get observations
		control.setError(new GaussianTimeError!dimensions(timeError, actualError, control.truth, params.integrator, &gen));
		control.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		//Do assimilation
		control.setLikelihood(new LikelihoodGetter!dimensions(control.observations, expectedError));
		control.getEnsembleTimeseries!false(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double controlRMSE = RMSE!dimensions(control.ensembleSeries, control.truth);
		
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood!dimensions(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter!dimensions(treatment.observations, expectedError);
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		//Deprecated; let the main loop handle file writing
		//File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", ", controlRMSE, ", ", treatmentRMSE);
		return [seed, observationInterval, timeError, controlRMSE, treatmentRMSE].to!string[1 .. $ - 1];
	} else if(params.config == RunConfigurations.CONTROL_RMSE) {
		Experiment!dimensions control = new Experiment!dimensions(params.integrator, params.controlAssimilator);
		//Get truth
		control.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		//Get observations
		control.setError(new GaussianTimeError!dimensions(timeError, actualError, control.truth, params.integrator, &gen));
		control.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		//Do assimilation
		control.setLikelihood(new LikelihoodGetter!dimensions(control.observations, expectedError));
		control.getEnsembleTimeseries!false(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double controlRMSE = RMSE!dimensions(control.ensembleSeries, control.truth);
		//Deprecated; let the main loop handle file writing
		//File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", ", controlRMSE);
		return [seed, observationInterval, timeError, controlRMSE].to!string[1 .. $ - 1];
	} else if(params.config == RunConfigurations.TREATMENT_RMSE) {
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood!dimensions(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter!dimensions(treatment.observations, expectedError);
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		//Deprecated; let the main loop handle file writing
		//File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", " , treatmentRMSE);
		return [seed, observationInterval, timeError, treatmentRMSE].to!string[1 .. $ - 1];
	} else if(params.config == RunConfigurations.INFERRED_TIME_ERROR) {
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood!dimensions(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter!dimensions(treatment.observations, expectedError);
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 2, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		DiscreteExperimentalLikelihood!dimensions treatmentLikelihood = cast(DiscreteExperimentalLikelihood!dimensions) treatment.likelihoodGetter;
		//Deprecated; let the main loop handle file writing
		//File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", ", treatmentLikelihood.expectedTime, ",", treatmentLikelihood.timeDeviation, ",,", treatmentLikelihood.timeLikelihood.to!string[1 .. $ - 1]);
		return ([seed, observationInterval, timeError, treatmentLikelihood.expectedTime, treatmentLikelihood.timeDeviation, 0] ~ treatmentLikelihood.timeLikelihood).to!string[1 .. $ - 1];
	} else {
		assert(0, "Invalid Run Configuration selected!");
	}
}

enum RunConfigurations: string {
	COMPARE_RMSE = "Compare RMSE",
	CONTROL_RMSE = "Control RMSE",
	TREATMENT_RMSE = "Treatment RMSE",
	INFERRED_TIME_ERROR = "Inferred Time Error"
}

void main() {
	RunConfigurations config = RunConfigurations.INFERRED_TIME_ERROR;
	string filename = "data/dataCollection/InferredTimeError13.csv";
	string logfile = "data/ExperimentLog.txt";
	string tag = "Parallelism test";
	StopWatch stopwatch  = StopWatch(AutoStart.no); //Used for benchmarking

	//We need to know these at compiletime in order to run the experiment in parallel
	const uint observationIntervalNumber = 1;
	const uint timeErrorNumber = 15;
	const uint errorNumber = 1;
	const uint trials = 1;

	bool logThisExperiment = true; //Set this to false if you don't want to write the experiment to the file
	double[] observationIntervals = [1];
	double[] timeErrors = [];
	foreach(i; 0..timeErrorNumber) {
		timeErrors ~= i * 0.001;
	}
	double[] errors = [0.05];
	//This will set up a number of random seeds
	//The first map statement will give different random seeds every program run
	//The second map statement will ensure that all program runs are the same
	//You can also set random seeds to those outputted by the program to replicate its results
	ulong[] seeds = iota(0, trials, 1)
					/*.map!(a => unpredictableSeed)*/
					.map!(a => cast(ulong)a)
					.array;
	//Package the parameters into one object
	Parameters!dimensions params = Parameters!dimensions(
		new Vector!(double, dimensions)(1), //The initial point of the truth
		0, //The initial time with which to associate the initial point
		40, //The time at which to stop the experiment
		0.01, //The length of each step of the integrator
		new RK4!dimensions(new Line!dimensions(new Vector!(double, 1)(25))), //The integrator used to return points from previous points, and its system
		0, //When to start observing
		40, //When to stop observing
		0, //When to create the ensemble
		40, //When to stop assimilating
		0.01, //The step for ensemble integration
		0.1,//The amount of time the ensemble is run before beginning to assimilate
		new Vector!(double, dimensions)(1), //The mean of the initial ensemble distribution
		new Vector!(double, dimensions)(0.1), //The standard deviation of the initial ensemble distribution
		20, //The size of the ensemble
		new EAKF!dimensions(), //The assimilation method for the control
		new EAKF!dimensions(), //The assimilation method for the treatment 
		-0.06, //The first time that is a valid time for observation relative to reported time
		0.06, //The last time that is a valid time for observation relative to reported time
		41, //The amount of different time intervals tested in experimental likelihood algorithm
		observationIntervals, //The intervals between observations that will be tested
		timeErrors, //The time error standard deviations that will be tested
		errors, //The observation error standard deviations that will be tested
		seeds, //The seeds for each of the trials
		filename, //The file to write the data to
		config, //The run configuration
		tag //User-specified experiment information
	);
	writeln(params);
	stopwatch.start();
	if(logThisExperiment) File(logfile, "a").writeln(params);
	if(config == RunConfigurations.COMPARE_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Time Error, Control RMSE, Treatment RMSE");
	} else if(config == RunConfigurations.CONTROL_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Time Error, Control RMSE");
	} else if(config == RunConfigurations.TREATMENT_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Time Error, Treatment RMSE");
	} else if(config == RunConfigurations.INFERRED_TIME_ERROR) {
		string binMiddles = iota(0, params.bins, 1).map!(a => params.minimumOffset + (params.maximumOffset - params.minimumOffset) / (params.bins * 2) + a * (params.maximumOffset - params.minimumOffset) / (params.bins)).array.to!string[1 .. $ - 1];
		File(filename, "a").writeln("Seed, Observation Interval, Time Error, Inferred Time Offset, Inferred Time Error, Inferred Time Offset Distribution: , ", binMiddles);
	}
	string[trials][errorNumber][timeErrorNumber][observationIntervalNumber] results;
	foreach(obsIntIndex, observationInterval; observationIntervals.parallel) {
		foreach(timeErrorIndex, timeError; timeErrors.parallel) {
			foreach(errorIndex, error; errors.parallel) {
				foreach(trialIndex, ref seed; seeds.parallel) {
					results[obsIntIndex][timeErrorIndex][errorIndex][trialIndex] =   
						run(params, observationInterval, timeError, new Vector!(double, dimensions)(error), Random(seed), seed);
				}
			}
		}
	}
	static foreach(obsIntIndex; 0..observationIntervalNumber) {
		static foreach(timeErrorIndex; 0..timeErrorNumber) {
			static foreach(errorIndex; 0..errorNumber) {
				static foreach(trialIndex; 0..trials) {
					File(filename, "a").writeln(results[obsIntIndex][timeErrorIndex][errorIndex][trialIndex]);
				}
			}
		}
	}
	stopwatch.stop();
	writeln("Finished ", 
			observationIntervalNumber * timeErrorNumber * errorNumber * trials, 
			" trials in ", stopwatch.peek.total!"msecs", " ms.");
	if(logThisExperiment) File(logfile, "a").writeln("Complete in " ~ stopwatch.peek.total!"msecs".to!string ~ " ms!\n");
}