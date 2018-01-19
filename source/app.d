import std.conv;
import std.stdio;
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
import systems.System;
import systems.Lorenz63;

void main() {
	//Declare experiment parameters
	//Universal
	const Vector startState = Vector(0, 0, 0); ///The initial point of the truth
	const double startTime = 0; ///The initial time with which to associate the initial point
	const double endTime = 200; ///The time at which to stop the experiment
	const double dt = 0.01; ///The length of each step of the integrator
	System system = new Lorenz63(); ///The dynamical system used as an environment for the experiment
	RK4 integrator = new RK4(system); ///The integrator used to return points from previous points
	//Getting observations
	Vector actualError = Vector(0.1, 0.1, 0.1); ///The standard deviation of the Gaussian error in space
	const double obsStartTime = 0; ///When to start observing
	const double obsEndTime = 0; ///When to stop observing
	const double timeError = 0.01; ///The standard deviation of the Gaussian error in time (may change to non-Gaussian in the future)
	const double observationInterval = 1; ///The time in between observations
	//Assimilation
	Vector expectedError = Vector(0.1, 0.1, 0.1); ///The a priori expected standard deviation for Gaussian space error
	const double ensembleStartTime = 0; ///When to create the ensemble
	const double ensembleEndTime = 200; ///When to stop assimilating
	const double ensembledt = 0.01; ///The step for ensemble integration
	const double spinup = 0; ///The amount of time the ensemble is run before beginning to assimilate
	const Vector ensembleGenesis = Vector(0, 0, 0); ///The mean of the initial ensemble distribution
	const Vector ensembleDeviation = Vector(1, 1, 1); ///The standard deviation of the initial ensemble distribution
	const int ensembleSize;
	Assimilator controlAssimilator = new EAKF(); ///The assimilation method for the control
	Assimilator experimentalAssimilator = new RHF(); ///The assimilation method for the treatment 
	const double minimumOffset = -0.2; ///The first time that is a valid time for observation relative to reported time
	const double maximumOffset = 0.2; ///The last time that is a valid time for observation relative to reported time
	const uint bins = 10; ///The amount of different time intervals tested in experimental likelihood algorithm
	//Experimental constants
	const string filename = "../data/testdata.csv"; ///The name of the file to write to

	writeln("Control:");
	Experiment control = new Experiment(integrator, controlAssimilator);
	control.getTruth(startState, startTime, endTime, dt);
	control.setError(new GaussianTimeError(timeError, actualError, control.truth, integrator));
	control.getObservations(obsStartTime, obsEndTime, observationInterval);
	control.setLikelihood(new LikelihoodGetter(control.observations, expectedError));
	control.getEnsembleTimeseries(
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
			treatment.observations, treatment.ensembleSeries, expectedError, integrator, minimumOffset, maximumOffset, bins
		)
	);
	treatment.getEnsembleTimeseries(
		ensembleStartTime, ensembleEndTime, ensembledt, spinup, new Ensemble(ensembleGenesis, ensembleSize, ensembleDeviation)
	);
	immutable double treatmentRMSE = RMSE(treatment.ensembleSeries, treatment.truth);
	writeln("Treatment RMSE is ", treatmentRMSE);
	write(filename, controlRMSE.to!string ~ ", " ~ treatmentRMSE.to!string);
}