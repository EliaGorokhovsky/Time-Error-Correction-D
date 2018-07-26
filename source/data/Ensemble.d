/**
 * The ensemble is essentially a Monte Carlo approximation of what we assume to be a Gaussian distribution
 * It represents our prior(forecast) as a set of potential system states
 * Ideally, we can assume that the true state and the ensemble members are from the same distribution
 */
module data.Ensemble;

import std.algorithm; //Used for map statements
import std.array; //Used to convert map statement outputs to arrays
import std.conv; //Used for casting to string
import std.range; //Used to generate stepped ranges of numbers
import mir.random; //Used to generate a random ensemble from a start point 
import mir.random.variable; //Used to generate a random ensemble from a start point
import utility.ArrayStats; //Used to find variance and standard deviation
import math.Vector; //Used to store state

/**
 * An ensemble of points (position vectors)
 * Functionally an array
 */
class Ensemble(uint dim) {

    Vector!(double, dim)[] members; ///The points of the ensemble

    /**
     * Returns a 2d array representing the transpose of the ensemble - arrays of its members' dimensional values
     */
    @property double[][dim] valueLists() {
        double[][dim] lists;
        foreach (component; this.members) {
            static foreach (i; 0..dim) {
                lists[i] ~= component[i];
            }
        }
        return lists;
    }

    /** 
     * Gets a Vector containing the mean values of the ensemble
     */
    @property Vector!(double, dim) eMean() {
        return new Vector!(double, dim)(this.valueLists[].map!(a => a.mean).array);
    }

    /**
     * Gets a vector containing variance for the ensemble members in each variable
     */
    @property Vector!(double, dim) eVariance() {
        return new Vector!(double, dim)(this.valueLists[].map!(a => a.variance!1).array);
    }

    /**
     * Gets a vector containing standard deviation for the ensemble members in each variable
     */
    @property Vector!(double, dim) eStandardDeviation() {
        return new Vector!(double, dim)(this.valueLists[].map!(a => a.standardDeviation!1).array);
    }

    /** 
     * Gets the size of the ensemble
     */
    @property ulong size() {
        return this.members.length;
    }

    /** 
     * Initializer for an ensemble
     * Generates ensemble with independent Gaussian variation from a base point 
     */
    this(Vector!(double, dim) base, int size, Vector!(double, dim) error, Random* gen) {
        foreach(i; 0..size) {
            double[dim] vars = iota(0, dim, 1).map!(a => NormalVariable!double(base[a], error[a])(*gen)).array;
            this.members ~= new Vector!(double, dim)(vars);
        }
    }

    /**
     * Constructor for an ensemble
     * Constructs an ensemble using a set of members
     */
    this(Vector!(double, dim)[] points) {
        this.members = points;
    }

    /**
     * Constructor for an ensemble
     * Constructs an ensemble using three lists of values that become xValues, yValues, and zValues
     */
    this(double[][dim] valueLists) {
        assert(valueLists[].all!(a => a.length == valueLists[0].length));
        double[dim] vectorComponents;
        foreach(i; 0..valueLists[0].length) {
            static foreach (j; 0..dim) {
                vectorComponents[j] = valueLists[j][i];
            }
            this.members ~= new Vector!(double, dim)(vectorComponents);
        }
    }

    /**
     * Adding a vector to an ensemble returns a new ensemble linearly shifted by that vector
     */
    Ensemble!dim opBinary(string op)(Vector!(double, dim) other) {
        static if (op == "+") return new Ensemble!dim(this.members.map!(a => a + other).array);
        else static if (op == "-") return new Ensemble!dim(this.members.map!(a => a - other).array);
        else return null;
    }

    /**
     * Scaling the ensemble should keep the mean the same but adjust the standard deviation
     */
    Ensemble!dim opBinary(string op)(double scalar) {
        static if (op == "*") return new Ensemble!dim(this.members.map!(a => this.eMean + (a - this.eMean) * scalar).array);
        else static if (op == "/") return new Ensemble!dim(this.members.map!(a => this.eMean + (a - this.eMean) / scalar).array);
        else return null;
    }

    /**
     * Handles shifting the ensemble without creating a new one
     */
    void opOpAssign(string op)(Vector!(double, dim) other) {
        static if (op == "+=") this.members = this.members.map!(a => a + other).array;
        else static if (op == "-=") this.members = this.members.map!(a => a - other).array;
    }

    /**
     * Handles scaling the ensemble without creating a new one
     */
    void opOpAssign(string op)(double scalar) {
        Vector!(double, dim) average = this.eMean;
        static if (op == "*=") this.members = this.members.map!(a => average + (a - average) * scalar).array;
        else static if (op == "/=") this.members = this.members.map!(a => average + (a - average) / scalar).array;
    }

    /**
     * Returns a string representation of the ensemble
     */
    override string toString() {
        return this.members.to!string;
    }

    /**
     * Returns a copy of the ensemble
     */
    Ensemble!dim copy() {
        return new Ensemble!dim(this.members.map!(a => new Vector!(double, dim)(a)).array);
    }

}

unittest {

    import std.stdio;
    
    writeln("\nUNITTEST: Ensemble");
    Vector!(double, 3) base = new Vector!(double, 3)(0);
    Random gen = Random(1);
    Ensemble!3 ensemble = new Ensemble!3(base, 20, new Vector!(double, 3)(0.01), &gen);
    writeln("Ensemble with base (0, 0, 0) and error (0.01, 0.01, 0.01) has x values ", ensemble.valueLists[0]);
    writeln("Ensemble mean in x is ", ensemble.eMean.x);
    writeln("Ensemble standard deviation in x is ", ensemble.eStandardDeviation.x);
    ensemble += new Vector!(double, 3)(2);
    writeln("Shifting the ensemble by 2 results in mean ", ensemble.eMean);
    ensemble *= 2;
    writeln("Scaling the ensemble by 2 results in standard deviation, ", ensemble.eStandardDeviation);

}