/**
 * A timeseries in dynamical systems research is any function of time
 * Here we store a list of states and map them to a list of times
 * We can use an integrator to interpolate and get the value of the function at an undefined time
 * We can also use linear interpolation to find it for vector timeseries
 * This is basically a data storage class
 */
module data.Timeseries;

import std.algorithm; //Used for sorting when finding dt
import std.math; //Used for approximate equality among doubles
import std.traits; //Used to check if the timeseries type is of a Vector or Ensemble template
import data.Ensemble; //Used as an input type for an overload of value() that interpolates
import integrators.Integrator; //Used to interpolate for getting value at undefined times
import math.Vector; //Used to define possible data storage objects that can be interpolated

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
     * Do not use this; favor value() instead
     */
    @property T[double] timeAssociate() {
        assert(this.members.length == this.times.length);
        T[double] association;
        //Fills the associative array
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
     * Does so by finding the time difference between the two smallest times in the timeseries
     * This should be true in most cases for this experiment, but be careful when using it elsewhere
     */
    @property double dt() {
        return this.times.sort[1] - this.times.sort[0];
    }

    /**
     * If the timeseries contains ensembles, returns a timeseries of their means
     */
    static if(__traits(isSame, TemplateOf!(T), Ensemble)) {
        @property Timeseries!(Vector!(double, dim)) meanSeries(uint dim)() {
            assert(this.members.length == this.times.length);
            Timeseries!(Vector!(double, dim)) means = new Timeseries!(Vector!(double, dim))();
            //Fils the new timeseries
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
    static if (__traits(isSame, TemplateOf!(T), Vector) || __traits(isSame, TemplateOf!(T), Ensemble)) {
        T value(uint dim)(double time, Integrator!dim integrator = null) {
            //If the time is defined, no need to use the integrator
            if(this.times.any!(a => a.approxEqual(time, 1e-6, 1e-6))) { 
                double newTime = this.times.find!(a => a.approxEqual(time, 1e-6, 1e-6))[0];
                return this.timeAssociate[newTime]; 
            }
            else {
                //Otherwise, find the last defined time before the given time
                ulong lastCountedTime = clamp(this.times.countUntil!"a > b"(time) - 1, 0, 100000);
                //And use a single integrator step to get to the given time from there
                double dt = time - this.times[cast(uint)lastCountedTime];
                return integrator(this.timeAssociate[this.times[cast(uint)lastCountedTime]], dt);
            }
        }
    }

    /**
     * Only works if the exact value exists in the times.
     */
    T valueAtTime(double time) {
        if(this.times.any!(a => a.approxEqual(time, 1e-6, 1e-6))) { 
            double newTime = this.times.find!(a => a.approxEqual(time, 1e-6, 1e-6))[0];
            return this.timeAssociate[newTime]; 
        }
        return this.members[0];
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

    writeln("\nUNITTEST: Timeseries");
    Timeseries!(Vector!(double, 3)) timeseries= new Timeseries!(Vector!(double, 3))();
    class Test : System!3 {
        override Vector!(double, 3) opCall(Vector!(double, 3) state) { return new Vector!(double, 3)(1); }
    }
    foreach(i; 0..10) {
        timeseries.add(i, new Vector!(double, 3)(i));
    }
    writeln("Times: ", timeseries.times);
    writeln("Time Series: ", timeseries.members);
    writeln("Time Series(2): ", timeseries.value(2));
    writeln("Time Series(2.5): ", timeseries.value(2.5, new RK4!3(new Test())));

}