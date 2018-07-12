module utility.ArrayStats;

import std.algorithm;
import std.array;
import std.math;
import math;

/**
 * Returns the mean of a dataset
 */
deprecated {
    double amean(double[] dataset) {
        return dataset.sum / dataset.length;
    }
}

/**
 * Returns the variance of a dataset, with ddof being the delta degrees of freedom
 * When ddof = 1 the measure is unbiased but converges more slowly
 */
double variance(int ddof = 0)(double[] dataset) {
    return reduce!((a, b) => a + pow(b - mean(dataset), 2) / (dataset.length - ddof))(0.0, dataset);
}

/**
 * Returns the standard deviation of a dataset, with ddof being the delta degrees of freedom
 * When ddof = 1 the measure is unbiased but converges more slowly
 */
double standardDeviation(int ddof = 0)(double[] dataset) {
    return sqrt(variance!ddof(dataset));
}

/**
 * Computes the covariance matrix of a set of data, where dim is the number of dimensions and ddof being delta degrees of freedom
 */
Matrix!(double, dim, dim) covariance(uint dim, uint datasetSize, int ddof = 0)(double[][] data) {
    Matrix!(double, dim, dim) covarianceMatrix = new Matrix!(double, dim, dim)();
    foreach (i; 0..dim) {
        data[i][] -= data[i][].mean;
    }
    foreach (i; 0..dim) {
        foreach(j; 0..dim) {
            //Dot product of the datasets gives covariance; covariance with itself is variance
            covarianceMatrix[i][j] = dot!(double, datasetSize)(new Vector!(double, datasetSize)(data[i]), new Vector!(double, datasetSize)(data[j])) / (datasetSize - ddof);
        }
    }
    return covarianceMatrix;
}

unittest {
    import std.stdio;

    writeln("\nUNITTEST: ArrayStats");
    writeln("Covariance matrix of [[1,2,3],[4,5,6],[7,8,9] is ", covariance!(3, 3, 1)([[1,2,3],[4,5,6],[7,8,9]]));
}