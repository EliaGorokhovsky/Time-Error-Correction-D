/**
 * Contains various storage structs that are used to pass experiment parameters
 */
module experiment.Parameters;

import std.conv; //Used to convert things to strings
import assimilation.Assimilator; //Used to store the assimilation methods that will be used
import integrators.Integrator; //Used to store how the system will be advanced in time
import math.Vector; //Used to store start states
import systems.System; //Used to store the system that will be assimilated
import app; //Used to know about run configurations

/**
 * A class (basically a struct) storing experimental parameters that define how this experiment will run
 * Contains information pertinent to the entire experiment
 */
struct Parameters(uint dim) {

    Vector!(double, dim) startState; ///The initial point of the truth
    double startTime; ///The initial time with which to associate the initial point
	double endTime; ///The time at which to stop the experiment
	double dt; ///The length of each step of the integrator
	Integrator!dim integrator; ///The integrator used to return points from previous points
	//Getting observations
	double obsStartTime; ///When to start observing
	double obsEndTime; ///When to stop observing
	//Assimilation
	double ensembleStartTime; ///When to create the ensemble
	double ensembleEndTime; ///When to stop assimilating
	double ensembledt; ///The step for ensemble integration
	double spinup; ///The amount of time the ensemble is run before beginning to assimilate
	Vector!(double, dim) ensembleGenesis; ///The mean of the initial ensemble distribution
	Vector!(double, dim) ensembleDeviation; ///The standard deviation of the initial ensemble distribution
	uint ensembleSize; //The size of the ensemble
	Assimilator!dim controlAssimilator; ///The assimilation method for the control
	Assimilator!dim experimentalAssimilator; ///The assimilation method for the treatment 
	double minimumOffset; ///The first time that is a valid time for observation relative to reported time
	double maximumOffset; ///The last time that is a valid time for observation relative to reported time
	uint bins; ///The amount of different time intervals tested in experimental likelihood algorithm
    //Trials
    double[] observationIntervals; ///The different observation intervals that were tested
    double[] timeErrors; ///The different time errors that were tested
    double[] errors; ///The different state errors that were tested
    ulong[] seeds; ///The different random seeds that were used for the trials
    string datafile; ///The name of the file that is written to
    RunConfigurations config; ///The type of experiment that is being run
    string tag; ///User-specified information about the experiment

    /**
     * Represents the parameters as a string in order to write them to a file
     */
    string toString() {
        return
        "Experiment:" ~
        "\nFilename: " ~ this.datafile ~
        "\nRun Configuration: " ~ this.config ~
        "\nTag: " ~ this.tag ~
        "\nDimensions: " ~ dim.to!string ~
        "\nObservation Intervals: " ~ this.observationIntervals.to!string ~
        "\nTime Errors: " ~ this.timeErrors.to!string ~
        "\nState Errors: " ~ this.errors.to!string ~
        "\nTrials: " ~ this.seeds.length.to!string ~
        "\nTotal Runs: " ~ (this.timeErrors.length * this.errors.length * this.observationIntervals.length * this.seeds.length).to!string ~
        "\nSeeds: " ~ this.seeds.to!string ~
        "\nTruth Start State: " ~ this.startState.toString ~ 
        "\nTruth Start Time: " ~ this.startTime.to!string ~
        "; Truth End Time: " ~ this.endTime.to!string ~
        "; Truth dt: " ~ this.dt.to!string ~ 
        "\nSystem: " ~ this.integrator.slope.toString ~
        "\nIntegrator: " ~ this.integrator.toString ~
        "\nObservations Start Time: " ~ this.obsStartTime.to!string ~ 
        "; Observations End Time: " ~ this.obsEndTime.to!string ~
        "\nEnsemble Start Time: " ~ this.ensembleStartTime.to!string ~
        "; Ensemble End Time: " ~ this.ensembleEndTime.to!string ~
        "; Ensemble dt: " ~ this.ensembledt.to!string ~ 
        "; Spin-up Time: " ~ this.spinup.to!string ~
        "\nEnsemble Creation Mean: " ~ this.ensembleGenesis.toString ~ 
        "; Ensemble Creation Standard Deviation: " ~ this.ensembleDeviation.toString ~ 
        "\nEnsemble Size: " ~ this.ensembleSize.to!string ~ 
        "\nControl Assimilator: " ~ this.controlAssimilator.toString ~ 
        "; Treatment Assimilator: " ~ this.experimentalAssimilator.toString ~
        "\nMinimum Time Offset: " ~ this.minimumOffset.to!string ~
        "; Maximum Time Offset: " ~ this.maximumOffset.to!string ~
        "; Time Histogram Bins: " ~ this.bins.to!string;
    }

}