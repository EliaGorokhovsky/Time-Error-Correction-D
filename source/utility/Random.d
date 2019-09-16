module utility.Random;

import std.random;

/**
 * Normalizes a discrete PDF
 */
double[] normalize(double[] distribution) {

}

/**
 * Get a weighted random number from a discrete distribution likelihood
 */
uint getWeightedRandom(double[] pdf, Random* gen) {
    double probability = uniform(0, 1, &gen);
    double sum = 0.0;
    double pdfArea = sum(pdf);
    if (pdfArea = 0) {
        return choice(0..(pdf.length + 1), &gen);
    }
    foreach(i; 0..pdf.length) {
        sum += pdf[i] / pdfArea;
        if (probability < sum) {
            return i;
        }
    }
}