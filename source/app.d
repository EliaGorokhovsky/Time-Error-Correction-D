import std.algorithm;
import std.array;
import std.conv;
import std.file;
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
import experiment.error.GaussianTimeError;
import integrators.RK4;
import systems.Circle;
import systems.System;
import systems.Lorenz63;

void run(double observationInterval, double timeError, Vector error) {
	//Declare experiment parameters
	//Universal
	Vector startState = Vector(1, 1, 1); ///The initial point of the truth
	const double startTime = 0; ///The initial time with which to associate the initial point
	const double endTime = 100; ///The time at which to stop the experiment
	const double dt = 0.01; ///The length of each step of the integrator
	Lorenz63 system = new Lorenz63(); ///The dynamical system used as an environment for the experiment
	RK4 integrator = new RK4(system); ///The integrator used to return points from previous points
	//Getting observations
	Vector actualError = error; ///The standard deviation of the Gaussian error in space
	const double obsStartTime = 0; ///When to start observing
	const double obsEndTime = 100; ///When to stop observing
	//Assimilation
	Vector expectedError = error; ///The a priori expected standard deviation for Gaussian space error
	const double ensembleStartTime = startTime; ///When to create the ensemble
	const double ensembleEndTime = endTime; ///When to stop assimilating
	const double ensembledt = dt; ///The step for ensemble integration
	const double spinup = 0.1; ///The amount of time the ensemble is run before beginning to assimilate
	const Vector ensembleGenesis = Vector(1, 1, 1); ///The mean of the initial ensemble distribution
	const Vector ensembleDeviation = Vector(1, 1, 1); ///The standard deviation of the initial ensemble distribution
	const int ensembleSize = 80;
	EAKF controlAssimilator = new EAKF(); ///The assimilation method for the control
	RHF experimentalAssimilator = new RHF(); ///The assimilation method for the treatment 
	const double minimumOffset = -0.1; ///The first time that is a valid time for observation relative to reported time
	const double maximumOffset = 0.1; ///The last time that is a valid time for observation relative to reported time
	const uint bins = 10; ///The amount of different time intervals tested in experimental likelihood algorithm
	//Experimental constants
	const string filename = "data/dataCollection/ContinuousTimeErrorData2.csv"; ///The name of the file to write to
	File file = File(filename, "a"); ///The file to be written to

	writeln("Control:");
	Experiment control = new Experiment(integrator, controlAssimilator);
	//Get truth
	control.getTruth(startState, startTime, endTime, dt);
	//Get observations
	control.setError(new GaussianTimeError(timeError, actualError, control.truth, integrator));
	control.getObservations(obsStartTime, obsEndTime, observationInterval);
	//Do assimilation
	control.setLikelihood(new LikelihoodGetter(control.observations, expectedError));
	control.getEnsembleTimeseries!false(
		ensembleStartTime, ensembleEndTime, ensembledt, spinup, new Ensemble(ensembleGenesis, ensembleSize, ensembleDeviation)
	);
	immutable double controlRMSE = RMSE(control.ensembleSeries, control.truth);
	writeln("Control RMSE is ", controlRMSE);
	writeln("Experiment:");
	Experiment treatment = new Experiment(integrator, experimentalAssimilator);
	treatment.getTruth(startState, startTime, endTime, dt);
	treatment.setError(new GaussianTimeError(timeError, actualError, treatment.truth, integrator));
	treatment.getObservations(obsStartTime, obsEndTime, observationInterval);
	treatment.setLikelihood(
		new DiscreteExperimentalLikelihood(
			treatment.observations, expectedError, integrator, minimumOffset, maximumOffset, bins
		)
	);
	treatment.getEnsembleTimeseries!true(
		ensembleStartTime, ensembleEndTime, ensembledt, spinup, new Ensemble(ensembleGenesis, ensembleSize, ensembleDeviation)
	);
	immutable double treatmentRMSE = RMSE(treatment.ensembleSeries, treatment.truth);
	writeln("Treatment RMSE is " ~ treatmentRMSE.to!string);
	file.writeln(observationInterval, ", ", timeError, ", ", error.x, ", ", controlRMSE, ", ", treatmentRMSE);
}

void main() {
	double[] observationIntervals = [1];
	double[] timeErrors = [];
	foreach(i; 0 .. 20) {
		timeErrors ~= cast(double)i / 200;
	}
	double[] errors = [0.1];
	uint trials = 5;
	File("data/dataCollection/ContinuousTimeErrorData2.csv", "a").writeln("Observation Interval, Time Error, State Error");
	//TODO: Make this clearer
	foreach(observationInterval; observationIntervals) {
		foreach(timeError; timeErrors) {
			foreach(error; errors) {
				foreach(i; 0..trials) {
					run(observationInterval, timeError, Vector(error, error, error));
				}
			}
		}
	}
}