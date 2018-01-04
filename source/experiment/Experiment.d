module experiment.Experiment;

import std.range;
import std.typecons;
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

    Timeseries!Vector truth;
    Timeseries!Vector observations;

    this(Integrator integrator) {
        this.integrator = integrator;
    }

    /**
     * Integrates a start point over the specified time interval
     * Ends at the last timestep before endTime
     */
    Timeseries!Vector getTruth(Vector point, double startTime, double endTime, double dt) {
        Timeseries!Vector truth = new Timeseries!Vector([point], [0]);
        foreach(i; iota(startTime, endTime, dt)) {
            point = integrator(point, dt);
            truth.add(i + dt, point);
        }
        this.truth = truth;
        return this.truth;
    }

    /**
     * Gets a set of observations 
     * Ends at the last interval before endTime
     * TODO
     */
    Timeseries!Vector getObservations(double startTime, double endTime, double interval, ErrorGenerator errorGen) {
        Timeseries!Vector observations = new Timeseries!Vector();
        foreach(i; iota(startTime, endTime, interval)) {
            observations.add(i, errorGen(i));
        }
        this.observations = observations;
        return this.observations;
    }

}

unittest {

    import std.stdio;
    import integrators.RK4;
    import experiment.error.GaussianError;

    writeln("\nUNITTEST: Experiment");
    class Test : System {
        override Vector opCall(Vector state) { return Vector(1, 1, 1); }
    }
    RK4 rk4 = new RK4(new Test());
    Experiment process = new Experiment(rk4);
    process.getTruth(Vector(0, 0, 0), 0, 10, 1);
    writeln("Integrating <1, 1, 1> from 0 to 10 returns ", process.truth.members);
    GaussianError error = new GaussianError(Vector(0.1, 0.1, 0.1), process.truth, rk4);
    process.getObservations(0, 10, 3, error);
    writeln("Observing every 3 seconds with std (0.1, 0.1, 0.1) returns ", process.observations.members);

}