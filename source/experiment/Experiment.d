module experiment.Experiment;

import std.range;
import std.typecons;
import data.Ensemble;
import data.Timeseries;
import data.Vector;
import experiment.Analytics;
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
        Timeseries!Vector truth = new Timeseries!Vector([point]);
        foreach(i; iota(startTime, endTime, dt)) {
            point = integrator(point, dt);
            truth.add(i, point);
        }
        this.truth = truth;
        return this.truth;
    }

    /**
     * Gets a set of observations 
     * Ends at the last interval before endTime
     * TODO
     */
    Timeseries!Vector getObservations(double interval, double startTime, double endTime) {
        return null;
    }

}

unittest {

    import std.stdio;
    import integrators.RK4;
    import systems.Lorenz63;

    writeln("\nUNITTEST: Process");
    RK4 rk4 = new RK4(new Lorenz63());
    Experiment experiment = new Experiment(rk4);
    writeln("Integrating <1, 1, 1> from 0 to 0.009 returns ", experiment.getTruth(Vector(1, 1, 1), 0, 0.01, 0.001));

}