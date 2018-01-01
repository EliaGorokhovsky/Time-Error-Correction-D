module assimilation.EAKF;

import std.algorithm;
import std.array;
import std.math;
import std.typecons;
import assimilation.Assimilator;
import data.Ensemble;
import data.Vector;
import utility.ArrayStats;
import utility.Regression;

/**
 * Ensemble Adjustment Kalman Filter as of Anderson(2003)
 */
class EAKF : Assimilator {

    Vector observation;
    Vector likelihood;

    /**
     * The constructor for an assimilator
     * Observation and likelihood are passed before assimilation
     */
    this(Vector observation, Vector likelihood) {
        this.observation = observation;
        this.likelihood = likelihood;
    }

    /**
     * Overloading to allow for calling the assimilator as a function
     * EAKF gets observation increments for a variable, regresses it onto all three variables, then repeats for the other variables
     */
    override Ensemble opCall(Ensemble prior) {
        Tuple!(Vector, Vector) posteriorMetrics = this.getPosterior(prior);
        Vector posteriorMean = posteriorMetrics[0];
        Vector posteriorSpread = posteriorMetrics[1];
        Ensemble output = prior.copy();
        double ySlope = regressionSlope(output.xValues, output.yValues);
        double zSlope = regressionSlope(output.xValues, output.zValues);
        double[] obsIncrements = this.getObservationIncrements(output.xValues, posteriorMean.x, posteriorSpread.x);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + ySlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + zSlope * obsIncrements[i]; }
        double xSlope = regressionSlope(output.yValues, output.xValues);
        zSlope = regressionSlope(output.yValues, output.zValues);
        obsIncrements = this.getObservationIncrements(output.yValues, posteriorMean.y, posteriorSpread.y);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + xSlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + zSlope * obsIncrements[i]; }
        xSlope = regressionSlope(output.zValues, output.xValues);
        ySlope = regressionSlope(output.zValues, output.yValues);
        obsIncrements = this.getObservationIncrements(output.zValues, posteriorMean.z, posteriorSpread.z);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + xSlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + ySlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        return output;
    }

    /**
     * Assuming normality, uses Bayes' Rule to find a posterior distribution given the instance's observation data.
     */
    Tuple!(Vector, Vector) getPosterior(Ensemble prior) {
        Vector standardDeviation = Vector(       //Posterior standard deviation is the square root of the inverse of the sum of the inverse squares of the prior and likelihood
            prior.eStandardDeviation.x == 0 || this.likelihood.x == 0 ? 
                0 : sqrt(1 / (pow(prior.eStandardDeviation.x, -2) + pow(this.likelihood.x, -2))),
            prior.eStandardDeviation.y == 0 || this.likelihood.y == 0 ? 
                0 : sqrt(1 / (pow(prior.eStandardDeviation.y, -2) + pow(this.likelihood.y, -2))),
            prior.eStandardDeviation.z == 0 || this.likelihood.z == 0 ? 
                0 : sqrt(1 / (pow(prior.eStandardDeviation.z, -2) + pow(this.likelihood.z, -2)))
        );
        Vector mean = Vector(
            prior.eStandardDeviation.x == 0 ? prior.eMean.x : this.likelihood.x == 0 ? this.observation.x :
                pow(standardDeviation.x, 2) * ((prior.eMean.x * pow(prior.eStandardDeviation.x, -2)) + (this.observation.x * pow(this.likelihood.x, -2))),
            prior.eStandardDeviation.y == 0 ? prior.eMean.y : this.likelihood.y == 0 ? this.observation.y :
                pow(standardDeviation.y, 2) * ((prior.eMean.y * pow(prior.eStandardDeviation.y, -2)) + (this.observation.y * pow(this.likelihood.y, -2))),
            prior.eStandardDeviation.z == 0 ? prior.eMean.z : this.likelihood.z == 0 ? this.observation.z :
                pow(standardDeviation.z, 2) * ((prior.eMean.z * pow(prior.eStandardDeviation.z, -2)) + (this.observation.z * pow(this.likelihood.z, -2)))
        );
        return tuple(mean, standardDeviation);
    }

    /**
     * Gets increments required to fit a set of values to a posterior
     */
    double[] getObservationIncrements(double[] priorValues, double posteriorMean, double posteriorSpread) {
        double meanDifference = posteriorMean - mean(priorValues);
        double spreadRatio;
        if (standardDeviation(priorValues) == posteriorSpread) {
            spreadRatio = 1;
        } else {
            spreadRatio = posteriorSpread / standardDeviation(priorValues);
        }
        return priorValues.map!(
            //   linear shift and stretch to equate prior distribution to posterior      
            a => posteriorMean + spreadRatio * (a - mean(priorValues)) - a
        ).array;
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: EAKF");
    EAKF eakf = new EAKF(Vector(0, 0, 0), Vector(1, 1, 1));
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
    writeln(eakf.getPosterior(test1));
    writeln(eakf.getPosterior(test2));
    writeln("Ensemble with mean ", test1.eMean, 
    " and standard deviation ", test1.eStandardDeviation, 
    " assimilated to <0, 0, 0> with likelihood standard deviation <1, 1, 1> returns an ensemble with mean ",
    result1.eMean, " and standard deviation ", result1.eStandardDeviation);
    writeln("Ensemble with mean ", test2.eMean, 
    " and standard deviation ", test2.eStandardDeviation, 
    " assimilated to <0, 0, 0> with likelihood standard deviation <1, 1, 1> returns an ensemble with mean ",
    result2.eMean, " and standard deviation ", result2.eStandardDeviation);

}