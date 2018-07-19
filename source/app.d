import std.algorithm;
import std.array;
import std.conv;
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
import data.Vector;
import data.Ensemble;
import experiment.Analytics;
import experiment.Experiment;
import experiment.Parameters;
import experiment.error.GaussianError;
import experiment.error.GaussianTimeError;
import experiment.error.UniformTimeError;
import integrators.Integrator;
import integrators.RK4;
import systems.TestSystem;
import systems.System;
import systems.Lorenz63;

void run(Parameters params, double observationInterval, double timeError, Vector error, Random gen, ulong seed) {
	
	Vector actualError = error; ///The standard deviation of the Gaussian error in space
	Vector expectedError = error; ///The a priori expected standard deviation for Gaussian space error
	if(params.config == RunConfigurations.COMPARE_RMSE) {
		writeln("Control:");
		Experiment control = new Experiment(params.integrator, params.controlAssimilator);
		//Get truth
		control.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		//Get observations
		control.setError(new GaussianTimeError(timeError, actualError, control.truth, params.integrator, &gen));
		control.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		//Do assimilation
		control.setLikelihood(new LikelihoodGetter(control.observations, expectedError));
		control.getEnsembleTimeseries!false(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double controlRMSE = RMSE(control.ensembleSeries, control.truth);
		writeln("Control RMSE for time error ", timeError, " is ", controlRMSE);

		writeln("Experiment:");
		Experiment treatment = new Experiment(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter(treatment.observations, expectedError);
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double treatmentRMSE = RMSE(treatment.ensembleSeries, treatment.truth);
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", ", controlRMSE, ", ", treatmentRMSE);
		writeln("Treatment RMSE for time error ", timeError, " is ", treatmentRMSE);
	} else if(params.config == RunConfigurations.CONTROL_RMSE) {
		writeln("Control:");
		Experiment control = new Experiment(params.integrator, params.controlAssimilator);
		//Get truth
		control.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		//Get observations
		control.setError(new GaussianTimeError(timeError, actualError, control.truth, params.integrator, &gen));
		control.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		//Do assimilation
		control.setLikelihood(new LikelihoodGetter(control.observations, expectedError));
		control.getEnsembleTimeseries!false(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double controlRMSE = RMSE(control.ensembleSeries, control.truth);
		writeln("Control RMSE for time error ", timeError, " is ", controlRMSE);
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", ", controlRMSE);
	} else if(params.config == RunConfigurations.TREATMENT_RMSE) {
		writeln("Experiment:");
		Experiment treatment = new Experiment(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter(treatment.observations, expectedError);
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		immutable double treatmentRMSE = RMSE(treatment.ensembleSeries, treatment.truth);
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", " , treatmentRMSE);
		writeln("Treatment RMSE for time error ", timeError, " is ", treatmentRMSE);
	} else if(params.config == RunConfigurations.INFERRED_TIME_ERROR) {
		writeln("Experiment:");
		Experiment treatment = new Experiment(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, observationInterval);
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood(
				treatment.observations, expectedError, params.integrator, params.minimumOffset, params.maximumOffset, params.bins, &gen
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter(treatment.observations, expectedError);
		treatment.getEnsembleTimeseries!true(
			params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		DiscreteExperimentalLikelihood treatmentLikelihood = cast(DiscreteExperimentalLikelihood) treatment.likelihoodGetter;
		File(params.datafile, "a").writeln(seed, ", ", observationInterval, ", ", timeError, ", ", treatmentLikelihood.expectedTime, ",", treatmentLikelihood.timeDeviation, ",,", treatmentLikelihood.timeLikelihood.to!string[1 .. $ - 1]);
		writeln("Inferred time error for time error ", timeError, " is ", treatmentLikelihood.timeDeviation);
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
	string filename = "data/dataCollection/InferredTimeError6.csv";
	string logfile = "data/ExperimentLog.txt";
	bool logThisExperiment = true; //Set this to false if you don't want to write the experiment to the file
	double[] observationIntervals = [0.5];
	double[] timeErrors = [];
	foreach(i; 0..50) {
		timeErrors ~= i * 0.001;
	}
	double[] errors = [0.1];
	//This will set up a number of random seeds
	//The first map statement will give different random seeds every program run
	//The second map statement will ensure that all program runs are the same
	//You can also set random seeds to those outputted by the program to replicate its results
	ulong[] seeds = iota(0, 1, 1)
					.map!(a => unpredictableSeed)
					/*.map!(a => cast(ulong)a)*/
					.array;
	//Package the parameters into one object
	Parameters params = Parameters(
		Vector(1, 1, 1), //The initial point of the truth
		0, //The initial time with which to associate the initial point
		20, //The time at which to stop the experiment
		0.01, //The length of each step of the integrator
		new TestSystem(), //The dynamical system used as an environment for the experiment
		new RK4(new TestSystem()), //The integrator used to return points from previous points
		0, //When to start observing
		20, //When to stop observing
		0, //When to create the ensemble
		20, //When to stop assimilating
		0.01, //The step for ensemble integration
		0.1,//The amount of time the ensemble is run before beginning to assimilate
		Vector(1, 1, 1), //The mean of the initial ensemble distribution
		Vector(0.0001, 0.0001, 0.0001), //The standard deviation of the initial ensemble distribution
		20, //The size of the ensemble
		new EAKF(), //The assimilation method for the control
		new EAKF(), //The assimilation method for the treatment 
		-0.1, //The first time that is a valid time for observation relative to reported time
		0.1, //The last time that is a valid time for observation relative to reported time
		20, //The amount of different time intervals tested in experimental likelihood algorithm
		observationIntervals, //The intervals between observations that will be tested
		timeErrors, //The time error standard deviations that will be tested
		errors, //The observation error standard deviations that will be tested
		seeds, //The seeds for each of the trials
		filename, //The file to write the data to
		config //The run configuration
	);
	writeln(params);
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
	foreach(observationInterval; observationIntervals) {
		foreach(timeError; timeErrors) {
			foreach(error; errors) {
				foreach(ref seed; seeds) {
					run(params, observationInterval, timeError, Vector(error, error, error), Random(seed), seed);
				}
			}
		}
	}
	if(logThisExperiment) File(logfile, "a").writeln("Complete!\n");
}