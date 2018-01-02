module utility.Normal;

import std.math;

/**
 * Uses a polynomial approximation of the inverse of the normal CDF to find the x for which the area under the normal function is equal to p
 */
double normInverse(double p) {
    //Originally by John Herrero.
    //http://home.online.no/~pjacklam/notes/invnorm
    double a1 = -39.69683028665376;
    double a2 =  220.9460984245205;
    double a3 = -275.9285104469687;
    double a4 =  138.357751867269;
    double a5 = -30.99479806614716;
    double a6 =  2.506628277459239;
    double b1 = -54.4760987982241;
    double b2 =  161.5858368580409;
    double b3 = -155.6989798598866;
    double b4 =  66.80131188771972;
    double b5 = -13.28068155288572;
    double c1 = -0.0778484002430293;
    double c2 = -0.3223964580411365;
    double c3 = -2.400758277161838;
    double c4 = -2.549732539343734;
    double c5 =  4.374664141464968;
    double c6 =  2.938163982698783;
    double d1 =  0.0084695709041462;
    double d2 =  0.3224671290700398;
    double d3 =  2.445134137142996;
    double d4 =  3.754408661907416;
    double pLow = 0.02425;
    double pHigh = 1 - pLow;
    double q;
    if(p < pLow) {
        q = sqrt(-2 * log(p));
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    } else if(p > pHigh) {
        q = sqrt(-2 * log(1 - p));
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    } else {
        q = p - 0.5;
        double r = q * q;
        return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
    }
}

/**
 * Find x such that the CDF of a normal distribution with certain mean and sd multiplied by alpha is p
 */
double weightedNormInverse(double alpha, double mean, double standardDeviation, double p) {
    double np = p / alpha;
    double x = normInverse(np);
    x = mean + x * standardDeviation;
    return x;
}