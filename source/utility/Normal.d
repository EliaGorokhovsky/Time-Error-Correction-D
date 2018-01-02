module utility.Normal;

import std.math;
import std.mathspecial;

/**
 * Find x such that the CDF of a normal distribution with certain mean and sd multiplied by alpha is p
 */
double weightedNormInverse(double alpha, double mean, double standardDeviation, double p) {
    double np = p / alpha;
    double x = normalDistributionInverse(np);
    x = mean + x * standardDeviation;
    return x;
}
