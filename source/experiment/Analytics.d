module experiment.Analytics;

import std.stdio;
import std.math;
import data.Ensemble;
import data.Timeseries;
import data.Vector;

/**
 * Finds the RMSE of a pair of timeseries
 * The RMSE is a measure of how accurate an ensemble is
 * RMSE is the square root of the mean of the squares deviations of the ensemble means
 */
double RMSE(Timeseries!Ensemble data, Timeseries!Vector truth) {
    if(data.members.length != truth.members.length) writeln("Data length: ", data.members.length, " Truth length: ", truth.members.length);
    assert(data.members.length == truth.members.length, "Unequal dataset lengths");
    double sumOfSquares = 0;
    foreach(i; 0..truth.members.length) {
        sumOfSquares += (truth.members[i].x - data.members[i].eMean.x).pow(2) + (truth.members[i].y - data.members[i].eMean.y).pow(2) + (truth.members[i].z - data.members[i].eMean.z).pow(2);
    }
    double squareMean = sumOfSquares / truth.members.length;
    return sqrt(squareMean);
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Analytics");
    Timeseries!Ensemble data = new Timeseries!Ensemble();
    Timeseries!Vector truth = new Timeseries!Vector();
    foreach(i; 0..10) {
        data.add(i, new Ensemble([Vector(3, 0, 0), Vector(0, 0, 0), Vector(0, 0, 0)]));
        truth.add(i, Vector(0, 0, 0));
    }
    writeln("Ensembles of ", data.members[0], " with Vector ", truth.members[0], " result in an RMSE of ", RMSE(data, truth));

}