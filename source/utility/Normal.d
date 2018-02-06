module utility.Normal;

import std.math;
import std.mathspecial;

/**
 * Find x such that the CDF of a normal distribution with certain mean and sd multiplied by alpha is p
 */
double weightedNormInverse(double alpha, double mean, double standardDeviation, double p) {
    double np = p / alpha;
    assert((1 - np) * (-np) < 0, "np is outside of bound");
    double x = normalDistributionInverse(np);
    assert(!isNaN(cast(float)x), "Normal(13): x is NaN");
    x = mean + x * standardDeviation;
    return x;
}

/**
 * Gets the probability of the occurrence of something under a value in a defined normal distribution
 * Actually the the area under the normal distribution from -inf. to val
 */
double normalIntegral(double val, double mu = 0, double sigma = 1) {
    return normalDistribution(val * sigma + mu);
}

/**
 * Gets the pprobability density of a value in a normal distribution
 * Actually the area under the normal distribution from val - x to val + x as x -> 0
 * This is calculated as:
 *           1                 -(x - mu)^2 / (2sigma^2)    
 *   -----------------   *   e
 *   sqrt(2pi*sigma^2)
 */
double normalVal(double val, double mu = 0, double sigma = 1) {
    return E.pow(- (val - mu).pow(2) / (2 * sigma * sigma)) / sqrt(2 * PI * sigma * sigma);
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Normal");
    writeln("Matplotlib returns normpdf:\n0,0,1 => 0.3989422804014327\n1,0,1 => 0.24197072451914337\n-1,1,2 => 0.12098536225957168");
    writeln("normalVal returns normalVal:\n0,0,1 => ", normalVal(0, 0, 1), "\n1,0,1 => ", normalVal(1, 0, 1), "\n-1,1,2 => ", normalVal(-1, 1, 2));

}
