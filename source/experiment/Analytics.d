module experiment.Analytics;

import std.algorithm;
import std.conv;
import std.math;
import std.range;
import std.stdio;
import std.traits;
import data.Ensemble;
import data.Timeseries;
import math.Vector;

/**
 * Finds the RMSE of a pair of timeseries
 * The RMSE is a measure of how accurate an ensemble is
 * RMSE is the square root of the mean of the squares deviations of the ensemble means
 */
double RMSE(uint dim)(Timeseries!(Ensemble!dim) data, Timeseries!(Vector!(double, dim)) truth) {
    assert(data.members.length == truth.members.length, "Unequal dataset lengths: " ~ data.members.length.to!string ~ " and " ~ truth.members.length.to!string);
    double sumOfSquares = 0;
    foreach(i; 0..truth.members.length) {
        sumOfSquares += iota(0, dim, 1).map!(a => (truth.members[i][a] - data.members[i].eMean[a]).pow(2)).sum;
    }
    double squareMean = sumOfSquares / truth.members.length;
    return sqrt(squareMean);
}

/**
 * Verifies if any elements in an array are NaN
 */
bool checkNaN(T)(T[] toCheck) {
    return toCheck.any!(a => a.isNaN);
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Analytics");
    Timeseries!(Ensemble!3) data = new Timeseries!(Ensemble!3)();
    Timeseries!(Vector!(double, 3)) truth = new Timeseries!(Vector!(double, 3))();
    foreach(i; 0..10) {
        data.add(i, new Ensemble([new Vector!(double, 3)([3, 0, 0]), new Vector!(double, 3)(0), new Vector!(double, 3)(0)]));
        truth.add(i, new Vector!(double, 3)(0));
    }
    writeln("Ensembles of ", data.members[0], " with Vector ", truth.members[0], " result in an RMSE of ", RMSE(data, truth));

}