module experiment.Analytics;

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
    assert(data.length == truth.length);
    double sumOfSquares = 0;
    foreach(i; 0..truth.length) {
        sumOfSquares += (truth[i].x - data[i].eMean.x).pow(2) + (truth[i].y - data[i].eMean.y).pow(2) + (truth[i].z - data[i].eMean.z).pow(2);
    }
    double squareMean = sumOfSquares / truth.length;
    return sqrt(squareMean);
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Analytics");
    Timeseries!Ensemble data = new Timeseries!Ensemble();
    Timeseries!Vector truth = new Timeseries!Vector();
    foreach(i; 0..10) {
        data.add(new Ensemble([Vector(3, 0, 0), Vector(0, 0, 0), Vector(0, 0, 0)]));
        truth.add(Vector(0, 0, 0));
    }
    writeln("Ensembles of ", data[0], " with Vector ", truth[0], " result in an RMSE of ", RMSE(data, truth));

}