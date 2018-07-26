/**
 * Contains the Ensemble Adjustment Kalman Filter, 
 * which was written by Dr. Jeff Anderson of NCAR
 * The EAKF finds a posterior by Gaussian multiplication
 * Then it adjusts the ensemble by shifting all the points so that their mean and spread
 * are equal to those of the posterior
 * I also account for covariance here using linear regression onto the other variables
 */
module assimilation.EAKF;

import std.algorithm; //Used for map functions
import std.array; //Used to convert map function outputs to arrays
import std.math; //Used to check approximate equality for doubles
import std.parallelism; //Used to quickly apply increments to ensemble members
import std.range; //Used to creat stepped forward ranges with iota
import std.typecons; //Used to handle tuples for posterior information outputs
import assimilation.Assimilator; //The parent class
import assimilation.likelihood.Likelihood; //Used for likelihood inputs
import data.Ensemble; //Used for the ensemble as inputs and outputs
import math.Vector; //Used for data storage 
import utility.ArrayStats; //Used for mean and standard deviation
import utility.Regression; //Used for regressing increments onto state variables

/**
 * Ensemble Adjustment Kalman Filter as of Anderson(2003)
 * Assimilates given Gaussian likelihood and assumed Gaussianity of prior
 * Adjusts prior ensemble with defined increments, rather than randomly like the Ensemble Kalman Filter (EnKF)
 */
class EAKF(uint dim) : Assimilator!dim {

    Vector!(double, dim) observation; ///The observation to assimilate
    Vector!(double, dim) likelihood; ///The standard deviation of a Gaussian likelihood with 0 covariance around the observation

    /**
     * The constructor for an assimilator
     * Observation and likelihood are passed before assimilation
     */
    this(Likelihood!dim likelihoodVars) {
        this.observation = likelihoodVars.gaussianMean;
        this.likelihood = likelihoodVars.gaussianDeviation;
    }

    /**
     * The empty constructor for an assimilator
     * We can construct an assimilator without any input, then set the likelihood later
     * if we intend to reuse it
     */
    this() {

    }

    /**
     * Sets the likelihood for the assimilator without creating a new one
     * Takes output of a LikelihoodGetter
     */
    override void setLikelihood(Likelihood!dim likelihoodVars) {
        this.observation = likelihoodVars.gaussianMean;
        this.likelihood = likelihoodVars.gaussianDeviation;
    }

    /**
     * Overloading to allow for calling the assimilator as a function
     * EAKF gets observation increments for a variable, regresses it onto all three variables, then repeats for the other variables
     */
    override Ensemble!dim opCall(Ensemble!dim prior) {
        //Gets the posterior in the form of a mean and standard deviation in 3d
        Tuple!(Vector!(double, dim), Vector!(double, dim)) posteriorMetrics = this.getPosterior(prior);
        Vector!(double, dim) posteriorMean = posteriorMetrics[0];
        Vector!(double, dim) posteriorSpread = posteriorMetrics[1];
        //Copies the prior ensemble to avoid changing it before doing calculations
        Ensemble!dim output = prior.copy();
        //Regress the observation increments for each variable onto the other variables
        double[dim] slopes;
        double[] obsIncrements;
        static foreach (i; 0..dim) {
            slopes = iota(0, dim, 1).map!(a => regressionSlope(output.valueLists[i], output.valueLists[a])).array;
            obsIncrements = this.getObservationIncrements(output.valueLists[i], posteriorMean[i], posteriorSpread[i]);
            static foreach (j; 0..dim) {
                foreach (index, increment; obsIncrements.parallel) {
                    output.members[index][j] = output.members[index][j] + slopes[j] * increment;
                }
            }
        }
        return output;
    }

    /**
     * Assuming normality, uses Bayes' Rule to find a posterior distribution given the instance's observation data.
     */
    Tuple!(Vector!(double, dim), Vector!(double, dim)) getPosterior(Ensemble!dim prior) {
        //Posterior standard deviation is the square root of the inverse of the sum of the inverse squares of the prior and likelihood
        //i.e. post = 1 / √((1/(prior^2)) + (1/(lik^2)))
        //in other words, the variance of the posterior is half the harmonic mean of the prior and likelihood variances
        //If either the likelihood or the prior standard deviation is zero then
        //the ensemble or the observation is 100% sure and then
        //the posterior should also have 0 variance
        double[dim] standardDeviations;
        static foreach (i; 0..dim) {
            standardDeviations[i] = prior.eStandardDeviation[i] == 0 || this.likelihood[i] == 0 ? 
                0 : sqrt(1 / (pow(prior.eStandardDeviation[i], -2) + pow(this.likelihood[i], -2)));
        }
        Vector!(double, dim) standardDeviation = new Vector!(double, dim)(standardDeviations);
        //The mean is the posterior variance times 
        //priorMean / priorStd^2 + obs / lik^2
        //There are several ugly ternary operators here:
        //If the prior standard deviation is zero then don't change it
        //If the likelihood standard deviation is zero just set the prior to it
        //If both are zero there's something wrong, so I don't change the ensemble
        double[dim] means;
        static foreach (i; 0..dim) {
            means[i] = prior.eStandardDeviation[i] == 0 ? prior.eMean[i] : this.likelihood[i] == 0 ? this.observation[i] :
                pow(standardDeviation[i], 2) * ((prior.eMean[i] * pow(prior.eStandardDeviation[i], -2)) + (this.observation[i] * pow(this.likelihood[i], -2)));
        }
        Vector!(double, dim) mean = new Vector!(double, dim)(means);
        return tuple(mean, standardDeviation);
    }

    /**
     * Gets increments required to fit a set of values to a posterior
     */
    double[] getObservationIncrements(double[] priorValues, double posteriorMean, double posteriorSpread) {
        double meanDifference = posteriorMean - mean(priorValues);
        double spreadRatio;
        //Get the ratio of the standard deviations of the ensemble members
        if (standardDeviation!1(priorValues) == posteriorSpread) {
            //This may look useless, but it only occurs when the prior values standard deviation is 0
            spreadRatio = 1;
        } else {
            spreadRatio = posteriorSpread / standardDeviation!1(priorValues);
        }
        //Adjust the ensemble
        return priorValues.map!(
            //linear shift and stretch to equate prior distribution to posterior      
            a => posteriorMean + spreadRatio * (a - mean(priorValues)) - a
        ).array;
    }

    /**
     * What to return when attempting to print this assimilator
     */
    override string toString() {
        return "EAKF";
    }

}

//Test several assimilations
unittest {

    import std.stdio;

    writeln("\nUNITTEST: EAKF");
    EAKF!3 eakf = new EAKF!3(new Likelihood(Vector!(double, 3)(0), Vector!(double, 3)(1)));
    //Standard deviation here should be zero, so ensemble should not change
    Ensemble test1 = new Ensemble(
        [1, 1, 1, 1, 1, 1 ,1, 1], 
        [1, 1, 1, 1, 1, 1 ,1, 1], 
        [1, 1, 1, 1, 1, 1 ,1, 1]
    );
    Ensemble test2 = new Ensemble(
        [0, 0.1, 0.2, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 1],
        [0, 0.1, 0.2, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 1],
        [0, 0.1, 0.2, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 1]
    );
    Ensemble result1 = eakf(test1);
    Ensemble result2 = eakf(test2);
    //writeln(eakf.getPosterior(test1));
    //writeln(eakf.getPosterior(test2));
    writeln("Ensemble with mean ", test1.eMean, 
    " and standard deviation ", test1.eStandardDeviation, 
    " assimilated to <0, 0, 0> with likelihood standard deviation <1, 1, 1> returns an ensemble with mean ",
    result1.eMean, " and standard deviation ", result1.eStandardDeviation);
    writeln("Ensemble with mean ", test2.eMean, 
    " and standard deviation ", test2.eStandardDeviation, 
    " assimilated to <0, 0, 0> with likelihood standard deviation <1, 1, 1> returns an ensemble with mean ",
    result2.eMean, " and standard deviation ", result2.eStandardDeviation);

}