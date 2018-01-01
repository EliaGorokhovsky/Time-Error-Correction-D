module utility.Regression;

import std.algorithm;
import std.array;
import std.math;
import dstats.regress;
import utility.ArrayStats;

/**
 * Calculates the slope of a regression line constructed from two same-size datasets
 */
double regressionSlope(double[] x, double[] y) {
    assert(x.length == y.length);
    int n = x.length;
    double[] xy;
    foreach (i; 0..n) {
        xy ~= x[i] * y[i];
    }
    double[] xSquared = x.map!(a => a * a).array;
    double[] ySquared = y.map!(a => a * a).array;
    double xDenom = n * xSquared.sum - x.sum * x.sum;
    double yDenom = n * ySquared.sum - y.sum * y.sum;
    if (xDenom == 0 || yDenom == 0) {
        return 1;
    }
    double r = (n * xy.sum - x.sum * y.sum) / sqrt(xDenom * yDenom);
    return r * standardDeviation(y) / standardDeviation(x);
}

/**
 * Calculates the intercept of a regression line constructed from two same-size datasets
 * Could be done more efficiently, but probably won't be used often
 */
double regressionIntercept(double[] x, double[] y) {
    return mean(y) - regressionSlope(x, y) * mean(x);
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Regression");
    double[] x = [9300, 10565, 15000, 15000, 17764, 57000, 65940, 73676, 77006, 93739, 146088, 153260];
    double[] y = [7100, 15500, 4400, 4400, 5900, 4600, 8800, 2000, 2750, 2550, 960, 1025];
    writeln("Regression for datasets:\n", x, "\n", y, "\nis y = ", regressionSlope(x, y), "x + ", regressionIntercept(x, y));

}