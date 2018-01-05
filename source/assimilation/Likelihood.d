module assimilation.Likelihood;

import std.typecons;
import data.Vector;

/**
 * A nullable data-storage class containing various ways to store a likelihood
 */ 
class Likelihood {

    Vector gaussianMean; ///A vector containing the mean of a gaussian likelihood
    Vector gaussianDeviation; ///A vector containing the standard deviation of a gaussian likelihood
    double[] xLikelihood; ///A list of probabilities in x dependent on the ensemble passed into a LikelihoodGetter
    double[] yLikelihood; ///A list of probabilities in y dependent on the ensemble passed into a LikelihoodGetter
    double[] zLikelihood; ///A list of probabilities in z dependent on the ensemble passed into a LikelihoodGetter

    /** 
     * A constructor for a likelihood
     * Used for Gaussian likelihoods with EAKF or similar
     */
    this(Vector gaussianMean, Vector gaussianDeviation) {
        this.gaussianMean = gaussianMean;
        this.gaussianDeviation = gaussianDeviation;
    }

    /**
     * A constructor for a likelihood
     * Used for non-Gaussian likelihoods with RHF
     */
    this(double[] xLikelihood, double[] yLikelihood, double[] zLikelihood) {
        this.xLikelihood = xLikelihood;
        this.yLikelihood = yLikelihood;
        this.zLikelihood = zLikelihood;
    }


}