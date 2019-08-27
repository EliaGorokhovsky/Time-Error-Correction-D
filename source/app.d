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
			params.testfile, params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
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
			params.testfile, params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		File(params.datafile, "a").writeln("'", seed, ", ", observationInterval, ",", error[0], ", ", timeError, ", ", controlRMSE, ", ", treatmentRMSE);
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
			params.testfile, params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 0, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double controlRMSE = RMSE!dimensions(control.ensembleSeries, control.truth);
		writeln("Control RMSE for time error ", timeError, " is ", controlRMSE);
		File(params.datafile, "a").writeln("'", seed, ", ", observationInterval, ",", error[0],  ", ", timeError, ", ", controlRMSE);
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
			params.testfile, params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 5, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		File(params.datafile, "a").writeln("'", seed, ", ", observationInterval, ",", error[0], ", ", timeError, ", " , treatmentRMSE, ", ");
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
			params.testfile, params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.spinup, 2, new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen)
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double RMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		DiscreteExperimentalLikelihood!dimensions treatmentLikelihood = cast(DiscreteExperimentalLikelihood!dimensions) treatment.likelihoodGetter;
		File(params.datafile, "a").writeln("'", seed, ", ", observationInterval, ",", error[0], ", ", timeError, ", ", treatmentLikelihood.expectedTime, ",", treatmentLikelihood.timeDeviation, ",", treatment.getTimeErrorRMSE(100), ",", treatment.getInferredTimeErrorRMSE(100), ",", treatment.getDirectionGuessRate(100), ",,", treatmentLikelihood.timeLikelihood.to!string[1 .. $ - 1]);
		writeln("Inferred time error for time error ", timeError, " is ", treatmentLikelihood.timeDeviation, " with observation time RMSE ", treatment.getTimeErrorRMSE(100), " and RMSE of calculated time error ", treatment.getInferredTimeErrorRMSE(100), " and direction guess rate ", treatment.getDirectionGuessRate(100), ". RMSE is ", RMSE);
	} else if (params.config == RunConfigurations.AVERAGE_DISTANCE) {
		writeln("Experiment:");
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		treatment.getTruth(params.startState, params.startTime, params.endTime, params.dt);
		treatment.setError(new GaussianTimeError!dimensions(timeError, actualError, treatment.truth, params.integrator, &gen));
		treatment.getObservations(params.obsStartTime, params.obsEndTime, params.obsEndTime - params.obsStartTime + 10);
		treatment.setLikelihood(new LikelihoodGetter!dimensions(treatment.observations, expectedError));
		treatment.getEnsembleTimeseries!false(
			params.testfile, params.ensembleStartTime, params.ensembleEndTime, params.ensembledt, params.ensembleEndTime, params.ensembleEndTime, new Ensemble!dimensions(params.ensembleGenesis, 5, params.ensembleDeviation, &gen)
		);
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		writeln("Treatment RMSE: ", treatmentRMSE);
	}
	else if (params.config == RunConfigurations.ALL) {
		Ensemble!dimensions ensemble = new Ensemble!dimensions(params.ensembleGenesis, params.ensembleSize, params.ensembleDeviation, &gen);

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
			params.testfile, 
			params.ensembleStartTime, 
			params.ensembleEndTime, 
			params.ensembledt, 
			params.spinup, 
			0,
			ensemble.copy() 
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double controlRMSE = RMSE!dimensions(control.ensembleSeries, control.truth);
		writeln("Control RMSE for time error ", timeError, " is ", controlRMSE);

		writeln("Experiment:");
		Experiment!dimensions treatment = new Experiment!dimensions(params.integrator, params.experimentalAssimilator);
		static if (verboseRun) writeln("Successfully initialized experiment.");
		treatment.setTruth(control.truth);
		static if (verboseRun) writeln("Successfully set true timeseries for treatment.");
		treatment.setObservations(control.observations, control.observationTimes);
		static if (verboseRun) writeln("Successfully set observations for treatment.");
		treatment.setLikelihood(
			new DiscreteExperimentalLikelihood!dimensions(
				treatment.observations, 
				expectedError, 
				params.integrator, 
				params.minimumOffset, 
				params.maximumOffset, 
				params.bins, 
				&gen, 
				//0, //known error
				//timeError //known error
			)
		);
		treatment.standardLikelihood = new LikelihoodGetter!dimensions(treatment.observations, expectedError);
		static if (verboseRun) writeln("Successfully set likelihood.");
		treatment.getEnsembleTimeseries!true(
			params.testfile, 
			params.ensembleStartTime, 
			params.ensembleEndTime, 
			params.ensembledt, 
			params.spinup, 
			5, 
			ensemble.copy()		
		);
		static if (verboseRun) writeln("Successfully ran ensemble.");
		immutable double treatmentRMSE = RMSE!dimensions(treatment.ensembleSeries, treatment.truth);
		DiscreteExperimentalLikelihood!dimensions treatmentLikelihood = cast(DiscreteExperimentalLikelihood!dimensions) treatment.likelihoodGetter;
		writeln("Treatment RMSE for time error ", timeError, " is ", treatmentRMSE);
		File(params.datafile, "a").writeln(
			"'", seed, ",",
			observationInterval, ",", 
			error[0], ",", 
			timeError, ",",
			treatmentLikelihood.expectedTime, ",", 
			treatmentLikelihood.timeDeviation, ",", 
			treatment.getTimeErrorRMSE(0), ",", 
			treatment.getInferredTimeErrorRMSE(0), ",", 
			treatment.getDirectionGuessRate(0), "," ,
			controlRMSE, ",",
			treatmentRMSE, ","
		);
	}

}

enum RunConfigurations: string {
	COMPARE_RMSE = "Compare RMSE",
	CONTROL_RMSE = "Control RMSE",
	TREATMENT_RMSE = "Treatment RMSE",
	INFERRED_TIME_ERROR = "Inferred Time Error",
	AVERAGE_DISTANCE = "Average Distance",
	ALL = "Run all tests"
}

void main() {
	RunConfigurations config = RunConfigurations.INFERRED_TIME_ERROR;
	string filename = "data/new-data/testing/onlySlope.csv";
	string testfilename = "data/new-data/testing/onlySlope.csv";
	string logfile = "data/new-data/testing/log-08-27.txt";
	string tag = "Test when slope is 20";
	StopWatch stopwatch = StopWatch(AutoStart.no);
	bool logThisExperiment = true; //Set this to false if you don't want to write the experiment to the logfile
	double[] observationIntervals = [1];
	double[] timeErrors = [0.05];
	/*foreach(i; 0..16) {
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
	ulong[] seeds = [unpredictableSeed];
					/*iota(0, 30, 1)
					.map!(a => unpredictableSeed)
					.map!(a => cast(ulong)a)
					.array;*/
					//[2184725850363998271uL, 13756863421133855013uL, 3923093028832091568uL, 9641058712962704131uL, 6750441739301156464uL, 12889439452371464117uL, 4309632604367651527uL, 18169754287967737752uL, 11881559119032824265uL, 6047906896178715925uL, 5463702583692681672uL, 15618503194663921397uL, 85746954789088986uL, 9711882800085160881uL, 9874236789977964086uL, 15837290847893381991uL, 8501143350149999419uL, 7021624231906446134uL, 1261576710274441696uL, 17612382575480320337uL];
					//[2628946988471595189uL, 15576251255597112670uL, 18254670456559473739uL, 13278908617683708147uL, 10571204722035012293uL, 11624264340839257134uL, 6739678969405065521uL, 4431202459603124066uL, 17728741792663068386uL, 15187366222404878511uL, 7244076727285745619uL, 9508980004864892384uL, 11462355250776930725uL, 15101795496735490436uL, 7060571836635634740uL, 15522236671386055511uL, 6883293110958960978uL, 9286824686537002161uL, 4526179969460858599uL, 12645961904936776150uL];
					/*[8370133083869323966uL]*/
					//[11533106597792568338, 4033616719091303501, 12266163737601295366, 16811336248043610353, 8624886549065854741, 6847397557744147839, 2882001111479385061, 2509445148698826592, 14131612931984432156, 16626537975641825998, 2295867116828041288, 11411048459360840025, 11040496366284795338, 5866271124652969060, 13704243046320118567, 6527747732338947434, 579817759488129967, 12866335293250994056, 4857661988290708509, 17364947456383381373, 2509457839729575406, 10604553522000053478, 15148233316620706376, 11309452247335074509, 12513937452241710350, 14174910347917962979, 17308290165121956148, 13140225654594366428, 17777330328431824263, 5417835852423128439]
					//[82029386284230530uL, 13963964892208834654uL, 9135740542362501819uL, 7031332123152652807uL, 15022117507543414108uL, 6622178055635864012uL, 1522060868933556120uL, 16100634437019494802uL, 8370133083869323966uL, 8296940596219331729uL];
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
		testfilename, //The file to write extra data to
		config, //The run configuration
		tag //User-specified experiment information
	);
	writeln(params);
	if(logThisExperiment) File(logfile, "a").writeln(params);
	stopwatch.start();
	if (config == RunConfigurations.COMPARE_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Control RMSE, Treatment RMSE");
	} else if (config == RunConfigurations.CONTROL_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Control RMSE");
	} else if (config == RunConfigurations.TREATMENT_RMSE) {
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Treatment RMSE");
	} else if (config == RunConfigurations.INFERRED_TIME_ERROR) {
		string binMiddles = iota(0, params.bins, 1).map!(a => params.minimumOffset + (params.maximumOffset - params.minimumOffset) / (params.bins * 2) + a * (params.maximumOffset - params.minimumOffset) / (params.bins)).array.to!string[1 .. $ - 1];
		File(filename, "a").writeln("Seed, Observation Interval, Observation Error, Time Error, Inferred Time Offset, Inferred Time Error, Observation Time RMSE, Adjusted Observation Time RMSE, Direction Guess Frequency, Inferred Time Offset Distribution: , ", binMiddles);
	} else if (config == RunConfigurations.AVERAGE_DISTANCE) {

	} else if (config == RunConfigurations.ALL) {
		File(filename, "a").writeln(
			"Seed, Observation Interval, Observation Error, Time Error, Inferred Time Offset, Inferred Time Error, Observation Time RMSE, Adjusted Observation Time RMSE, Direction Guess Rate, Control RMSE, Treatment RMSE"
		);
	}
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