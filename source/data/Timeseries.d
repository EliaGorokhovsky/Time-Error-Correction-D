module data.Timeseries;

import std.algorithm;
import std.math;
import data.Ensemble;
import data.Vector;
import integrators.Integrator;

/**
 * Essentially an array of something
 * Describes the state of an object over a number of steps
 * Contains a variety of useful analysis methods
 */
class Timeseries(T) {

    T[] members; ///The list of the object's states 
    double[] times; ///The list of times corresponding to the object's states

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
     * Gets the amount of members in the timeseries
     */
    @property ulong length() {
        return this.members.length;
    }

    /**
     * Returns a time increment, assuming difference between successive entries is constant
     * This should be true in most cases for this experiment, but be careful when using it elsewhere
     */
    @property double dt() {
        return this.times.sort[1] - this.times.sort[0];
    }

    /**
     * If the timeseries contains ensembles, returns a timeseries of their means
     */
    static if(is(T == Ensemble)) {
        @property Timeseries!Vector meanSeries() {
            assert(this.members.length == this.times.length);
            Timeseries!Vector means = new Timeseries!Vector();
            foreach(i; 0..this.members.length) {
                means.add(this.times[i], this.members[i].eMean);
            }
            return means;
        }
    }

    /**
     * Constructs a timeseries from an existing list of states
     */
    this(T[] members, double[] times) {
        this.members = members;
        this.times = times;
    }

    /**
     * Constructs a timeseries with empty arrays for time and state
     */
    this() {

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
        if(this.times.any!(a => a.approxEqual(time, 1e-6, 1e-6)))) { 
            double newTime = this.times.find!(a => a.approxEqual(time, 1e-6, 1e-6)))[0];
            return this.timeAssociate[newTime]; 
        }
        else {
            ulong lastCountedTime = (this.times.countUntil!"a > b"(time) > 0)? this.times.countUntil!"a > b"(time) - 1 : 0;
            double dt = time - this.times[cast(uint)lastCountedTime];
            return integrator(this.timeAssociate[this.times[cast(uint)lastCountedTime]], dt);
        }
    }

    /**
     * Gets the nth value of the timeseries
     */
    T value(uint index) {
        return this.members[index];
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
    writeln("Time Series(2): ", timeseries.value(2));
    writeln("Time Series(2.5): ", timeseries.value(2.5, new RK4(new Test())));

}