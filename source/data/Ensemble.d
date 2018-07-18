/**
 * The ensemble is essentially a Monte Carlo approximation of what we assume to be a Gaussian distribution
 * It represents our prior(forecast) as a set of potential system states
 * Ideally, we can assume that the true state and the ensemble members are from the same distribution
 */
module data.Ensemble;

import std.algorithm; //Used for map statements
import std.array; //Used to convert map statement outputs to arrays
import std.conv; //Used for casting to string
import mir.random; //Used to generate a random ensemble from a start point 
import mir.random.variable; //Used to generate a random ensemble from a start point
import data.Vector; //Used to represent an ensemble member
import utility.ArrayStats; //Used to find variance and standard deviation

/**
 * An ensemble of points (position vectors)
 * Functionally an array
 */
class Ensemble {

    Vector[] members; ///The points of the ensemble

    /**
     * Gets an array of all of the x-values of the points the ensemble
     */
    @property double[] xValues() {
        return members.map!(a => a.x).array;
    }

    /**
     * Sets the x-values independently of the other variables
     */
    @property void xValues(double[] input) {
        assert(input.length == this.members.length);
        foreach(i; 0..input.length) {
            this.members[i].x = input[i];
        }
    }

    /**
     * Gets an array of all of the y-values of the points of the ensemble
     */
    @property double[] yValues() {
        return members.map!(a => a.y).array;
    }

    /**
     * Sets the y-values independently of the other variables
     */
    @property void yValues(double[] input) {
        assert(input.length == this.members.length);
        foreach(i; 0..input.length) {
            this.members[i].y = input[i];
        }
    }

    /**
     * Gets an array of all of the z-values of the points of the ensemble
     */
    @property double[] zValues() {
        return members.map!(a => a.z).array;
    }

    /**
     * Sets the z-values independently of the other variables
     */
    @property void zValues(double[] input) {
        assert(input.length == this.members.length);
        foreach(i; 0..input.length) {
            this.members[i].z = input[i];
        }
    }

    /** 
     * Gets a Vector containing the mean values of the ensemble
     */
    @property Vector eMean() {
        return Vector(
            mean(this.xValues), 
            mean(this.yValues), 
            mean(this.zValues)
        );
    }

    /**
     * Gets a vector containing variance for the ensemble members in each variable
     */
    @property Vector eVariance() {
        return Vector(
            variance!1(this.xValues),
            variance!1(this.yValues),
            variance!1(this.zValues)
        );
    }

    /**
     * Gets a vector containing standard deviation for the ensemble members in each variable
     */
    @property Vector eStandardDeviation() {
        return Vector(
            standardDeviation!1(this.xValues),
            standardDeviation!1(this.yValues),
            standardDeviation!1(this.zValues)
        );
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
    this(Vector base, int size, Vector error, Random* gen) {
        auto normalX = NormalVariable!double(base.x, error.x);
        auto normalY = NormalVariable!double(base.y, error.y);
        auto normalZ = NormalVariable!double(base.z, error.z);
        foreach(i; 0..size) {
            this.members ~= Vector(normalX(*gen), normalY(*gen), normalZ(*gen));
        }
    }

    /**
     * Constructor for an ensemble
     * Constructs an ensemble using a set of members
     */
    this(Vector[] points) {
        this.members = points;
    }

    /**
     * Constructor for an ensemble
     * Constructs an ensemble using three lists of values that become xValues, yValues, and zValues
     */
    this(double[] xValues, double[] yValues, double[] zValues) {
        assert(xValues.length == yValues.length && yValues.length == zValues.length);
        foreach(i; 0..xValues.length) {
            this.members ~= Vector(xValues[i], yValues[i], zValues[i]);
        }
    }

    /**
     * Adding a vector to an ensemble returns a new ensemble linearly shifted by that vector
     */
    Ensemble opBinary(string op)(Vector other) {
        static if (op == "+") return new Ensemble(this.members.map!(a => a + other).array);
        else static if (op == "-") return new Ensemble(this.members.map!(a => a - other).array);
        else return null;
    }

    /**
     * Scaling the ensemble should keep the mean the same but adjust the standard deviation
     */
    Ensemble opBinary(string op)(double scalar) {
        static if (op == "*") return new Ensemble(this.members.map!(a => this.eMean + (a - this.eMean) * scalar).array);
        else static if (op == "/") return new Ensemble(this.members.map!(a => this.eMean + (a - this.eMean) / scalar).array);
        else return null;
    }

    /**
     * Handles shifting the ensemble without creating a new one
     */
    void opOpAssign(string op)(Vector other) {
        static if (op == "+=") this.members = this.members.map!(a => a + other).array;
        else static if (op == "-=") this.members = this.members.map!(a => a - other).array;
    }

    /**
     * Handles scaling the ensemble without creating a new one
     */
    void opOpAssign(string op)(double scalar) {
        Vector average = this.eMean;
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
    Ensemble copy() {
        return new Ensemble(this.members.dup);
    }

}

unittest {

    import std.stdio;
    
    writeln("\nUNITTEST: Ensemble");
    Vector base = Vector(0, 0, 0);
    Random gen = Random(unpredictableSeed);
    Ensemble ensemble = new Ensemble(base, 20, Vector(0.01, 0.01, 0.01), &gen);
    writeln("Ensemble with base (0, 0, 0) and error (0.01, 0.01, 0.01) has x values ", ensemble.xValues);
    writeln("Ensemble mean in x is ", ensemble.eMean.x);
    writeln("Ensemble standard deviation in x is ", ensemble.eStandardDeviation.x);
    ensemble += Vector(2, 2, 2);
    writeln("Shifting the ensemble by 2 results in mean ", ensemble.eMean);
    ensemble *= 2;
    writeln("Scaling the ensemble by 2 results in standard deviation, ", ensemble.eStandardDeviation);

}