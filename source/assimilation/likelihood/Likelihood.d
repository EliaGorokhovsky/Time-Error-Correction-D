/**
 * A storage class for assimilation likelihoods
 * This is necessary so that I can take in both discrete and continuous likelihoods on demand
 * It is a class because D structs are not nullable - it is advantageous for likelihood to be able to be null
 * as a placeholder
 * A likelihood is the probability distribution of a true state given an observation
 */
module assimilation.likelihood.Likelihood;

import math.Vector; //Used for data storage

/**
 * A nullable data-storage class containing various ways to store a likelihood
 */ 
class Likelihood(uint dim) {

    Vector!(double, dim) gaussianMean; ///A vector containing the mean of a gaussian likelihood
    Vector!(double, dim) gaussianDeviation; ///A vector containing the standard deviation of a gaussian likelihood
    double[dim][] likelihoods; ///A set of arrays corresponding to a dimension of the model (number of observed variables)

    /** 
     * A constructor for a likelihood
     * Used for Gaussian likelihoods with EAKF or similar
     */
    this(Vector!(double, dim) gaussianMean, Vector!(double, dim) gaussianDeviation) {
        this.gaussianMean = gaussianMean;
        this.gaussianDeviation = gaussianDeviation;
    }

    /**
     * A constructor for a likelihood
     * Used for non-Gaussian likelihoods with RHF
     */
    this(double[dim][] likelihoods) {
        this.likelihoods = likelihoods;
    }


}