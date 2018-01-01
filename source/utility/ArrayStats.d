module utility.ArrayStats;

import std.algorithm;
import std.array;
import std.math;

/**
 * Returns the mean of a dataset
 */
double mean(double[] dataset) {
    return dataset.sum / dataset.length;
}

/**
 * Returns the variance of a dataset, with ddof being the delta degrees of freedom
 */
double variance(double[] dataset, int ddof = 0) {
    return reduce!((a, b) => a + pow(b - mean(dataset), 2) / (dataset.length - ddof))(0.0, dataset);
}

/**
 * Returns the standard deviation of a dataset, with ddof being the delta degrees of freedom
 */
double standardDeviation(double[] dataset, int ddof = 0) {
    return sqrt(variance(dataset, ddof));
}