module data.Ensemble;

import std.algorithm;
import std.array;
import std.conv;
import std.math;
import mir.random;
import mir.random.variable;
import data.Vector;

/**
 * An ensemble of points (position vectors)
 * Functionally an array
 */
class Ensemble {

    Vector[] members;

    /**
     * Gets an array of all of the x-values of the points the ensemble
     */
    @property double[] xValues() {
        return members.map!(a => a.x).array;
    }

    /**
     * Gets an array of all of the y-values of the points of the ensemble
     */
    @property double[] yValues() {
        return members.map!(a => a.y).array;
    }

    /**
     * Gets an array of all of the z-values of the points of the ensemble
     */
    @property double[] zValues() {
        return members.map!(a => a.y).array;
    }

    /** 
     * Gets a Vector containing the mean values of the ensemble
     */
    @property Vector mean() {
        return Vector(
            this.xValues.sum/this.xValues.length, 
            this.yValues.sum/this.yValues.length, 
            this.zValues.sum/this.zValues.length
        );
    }

    /**
     * Gets a vector containing variance for the ensemble members in each variable
     */
    @property Vector variance() {
        return Vector(
            reduce!((a, b) => a + pow(b - this.mean.x, 2) / this.xValues.length)(0.0, this.xValues),
            reduce!((a, b) => a + pow(b - this.mean.y, 2) / this.yValues.length)(0.0, this.yValues),
            reduce!((a, b) => a + pow(b - this.mean.z, 2) / this.zValues.length)(0.0, this.zValues)
        );
    }

    /**
     * Gets a vector containing standard deviation for the ensemble members in each variable
     */
    @property Vector standardDeviation() {
        return Vector(sqrt(this.variance.x), sqrt(this.variance.y), sqrt(this.variance.z));
    }

    /** 
     * Initializer for an ensemble
     * Generates ensemble with independent Gaussian variation from a base point 
     */
    this(Vector base, int size, Vector error) {
        auto gen = Random(unpredictableSeed);
        auto normalX = NormalVariable!double(base.x, error.x);
        auto normalY = NormalVariable!double(base.y, error.y);
        auto normalZ = NormalVariable!double(base.z, error.z);
        foreach(i; 0..size) {
            this.members ~= Vector(normalX(gen), normalY(gen), normalZ(gen));
        }
    }

    /**
     * Initializer for an ensemble
     * Constructs an ensemble using a set of members
     */
    this(Vector[] points) {
        this.members = points;
    }

}

unittest {

    import std.stdio;
    
    writeln("\nUNITTEST: Ensemble");
    Vector base = Vector(0, 0, 0);
    Ensemble ensemble = new Ensemble(base, 20, Vector(0.01, 0.01, 0.01));
    writeln("Ensemble with base (0, 0, 0) and error (0.01, 0.01, 0.01) has x values ", ensemble.xValues);
    writeln("Ensemble mean in x is ", ensemble.mean.x);
    writeln("Ensemble standard deviation in x is ", ensemble.standardDeviation.x);

}