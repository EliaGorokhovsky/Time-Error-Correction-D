module data.Ensemble;

import std.algorithm;
import std.array;
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

    @property Vector mean() {
        return Vector(this.xValues.sum/this.xValues.length, this.yValues.sum/this.yValues.length, this.zValues.sum/this.zValues.length);
    }

    /** 
     * Initializer for an ensemble
     * Generates ensemble with independent Gaussian variation from a base point 
     * TODO
     */
    this(Vector base, int size, Vector error) {

    }

    /**
     * Initializer for an ensemble
     * Constructs an ensemble using a set of members
     */
    this(Vector[] points) {
        this.members = points;
    }

}