module data.Timeseries;

import std.algorithm;
import integrators.Integrator;

/**
 * Essentially an array of something
 * Describes the state of an object over a number of steps
 * Contains a variety of useful analysis methods
 */
class Timeseries(T) {

    T[] members; ///The list of the object's states 
    double[] times; ///The list of times corresponding to the object's states

    alias members this;

    /**
     * Returns an associative array associating time with state
     */
    @property T[double] timeAssociate() {
        assert(this.members.length == this.times.length);
        T[double] association;
        foreach(i; 0..this.members.length) {
            association[this.times[i]] = this.members[i];
        }
        return association;
    }

    /**
     * Constructs a timeseries from an existing list of states
     */
    this(T[] members = null) {
        this.members = members;
    }

    /**
     * Adds an element to the timeseries
     */
    void add(double time, T state) {
        this.times ~= time;
        this.members ~= state;
    }

    /**
     * Assuming times are sorted, finds the value of the state at a given time
     * If there is no corresponding state in the timeseries, use an integrator
     * TODO: Interpolate instead
     */
    T value(double time, Integrator integrator = null) {
        if(this.times.canFind(time)) { 
            return this.timeAssociate[time]; 
        }
        else {
            int lastCountedTime = this.times.countUntil!"a > b"(time) - 1;
            double dt = time - this.times[lastCountedTime];
            return integrator(this.timeAssociate[this.times[lastCountedTime]], dt);
        }
    }

}

unittest {

    import std.stdio;
    import systems.System;
    import integrators.RK4;
    import data.Vector;

    writeln("\nUNITTEST: Timeseries");
    Timeseries!Vector timeseries= new Timeseries!Vector();
    class Test : System {
        override Vector opCall(Vector state) { return Vector(1, 1, 1); }
    }
    foreach(i; 0..10) {
        timeseries.add(i, Vector(i, i, i));
    }
    writeln("Times: ", timeseries.times);
    writeln("Time Series: ", timeseries.members);
    writeln("Time Series(2): ", timeseries.value(2, new RK4(new Test())));
    writeln("Time Series(2.5): ", timeseries.value(2.5, new RK4(new Test())));

}