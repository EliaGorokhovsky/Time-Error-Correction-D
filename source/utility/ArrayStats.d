module utility.ArrayStats;

import std.algorithm;
import std.array;
import std.math;
import std.parallelism;
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
Matrix!(double, dim, dim) covariance(uint dim, int ddof = 0)(double[][] data) {
    Matrix!(double, dim, dim) covarianceMatrix = new Matrix!(double, dim, dim)();
    foreach (i; 0..dim) {
        data[i][] -= data[i][].mean;
    }
    foreach (i; 0..dim) {
        foreach(j; 0..dim) {
            assert(data[i].length == data[j].length, "Covariance: Dataset lengths unequal");
            //Dot product of the datasets gives covariance; covariance with itself is variance
            double productSum = 0;
            foreach (k; 0..data[i].length) {
                productSum += data[i][k] * data[j][k];
            }
            covarianceMatrix[i][j] = productSum / (data[i].length - ddof);
        }
    }
    return covarianceMatrix;
}

/**
 * Gets the mean of a discretely defined pdf represented as a list
 * Support is the domain of the function
 */
double PDFMean(double[] support, double[] probabilities) {
    assert(support.length == probabilities.length, "Function does not match its support");
    double[] supportCopy = support.dup;
    foreach(index, ref component; supportCopy.parallel) {
        component *= probabilities[index];
    }
    return supportCopy.sum / probabilities.sum;
}

/**
 * Gets the variance of a discretely defined pdf represented as a list
 * Support is the domain of the function
 */
double PDFVariance(double[] support, double[] probabilities) {
    assert(support.length == probabilities.length, "Function does not match its support");
    double[] supportCopy = support.dup;
    double mean = PDFMean(support, probabilities);
    foreach(index, ref component; supportCopy.parallel) {
        component = probabilities[index] * (component - mean).pow(2);
    }
    return supportCopy.sum / probabilities.sum;
}



unittest {
    
    import std.stdio;

    writeln("\nUNITTEST: ArrayStats");
    writeln("Covariance matrix of [[1,2,3],[4,5,6],[7,8,9] is ", covariance!(3, 1)([[1,2,3],[4,5,6],[7,8,9]]));
}