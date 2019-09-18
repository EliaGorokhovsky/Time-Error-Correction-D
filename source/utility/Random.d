module utility.Random;

import std.algorithm;
import std.conv;
import std.range;
import mir.random;
import mir.random.variable;

/**
 * Get a weighted random number from a discrete distribution likelihood
 */
uint getWeightedRandom(double[] pdf, Random* gen) {
    assert(gen !is null);
    double pdfArea = pdf.sum();
    if (pdfArea == 0) {
        auto uniformVar = UniformVariable!int(0, (pdf.length - 1).to!int);
        return uniformVar(*gen);
    }
    auto uniformVar = UniformVariable!double(0, 1);
    double probability = uniformVar(*gen);
    double cumulativeProbability = 0.0;
    foreach(i; 0..pdf.length) {
        cumulativeProbability += pdf[i] / pdfArea;
        if (probability < cumulativeProbability) {
            return i.to!uint;
        }
    }
    assert(0);
}