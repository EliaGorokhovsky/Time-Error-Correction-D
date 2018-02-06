module assimilation.RHF;

import std.algorithm;
import std.array;
import std.mathspecial;
import std.range;
import std.stdio;
import std.typecons;
import assimilation.Assimilator;
import assimilation.likelihood.Likelihood;
import data.Ensemble;
import data.Vector;
import utility.ArrayStats;
import utility.Math;
import utility.NDVector;
import utility.Normal;
import utility.Regression;
import utility.Sort;

/** 
 * The Rank Histogram Filter (Anderson 2010) to apply ensemble filtering with discretely defined likelihood distributions
 */
class RHF : Assimilator {

    double[] xLikelihood; ///The probabilities for the x-values of the ensemble
    double[] yLikelihood; ///The probabilities for the y-values of the ensemble
    double[] zLikelihood; ///The probabilities for the z-values of the ensemble
    bool rectangularQuadrature; ///Whether rectangular or trapezoidal quadrature will be used to calculate probabilities with likelihood

    this(Likelihood likelihood, bool rectangularQuadrature = true) {
        this.xLikelihood = likelihood.xLikelihood;
        this.yLikelihood = likelihood.yLikelihood;
        this.zLikelihood = likelihood.zLikelihood;
        this.rectangularQuadrature = rectangularQuadrature;

    }

    this() {

    }

    /**
     * Sets the likelihood for the assimilator without making a new one
     * Takes the output of a likelihoodGetter
     */
    override void setLikelihood(Likelihood likelihood) {
        this.xLikelihood = likelihood.xLikelihood;
        this.yLikelihood = likelihood.yLikelihood;
        this.zLikelihood = likelihood.zLikelihood;        
    }

    /**
     * Overloading to allow for calling the assimilator as a function
     * This RHF gets observation increments for a variable, regresses it onto all three variables, then repeats for the other variables
     */
    override Ensemble opCall(Ensemble prior) {
        Ensemble output = prior.copy();
        double ySlope = regressionSlope(output.xValues, output.yValues);
        double zSlope = regressionSlope(output.xValues, output.zValues);
        double[] obsIncrements = this.getObservationIncrements(output.xValues, this.xLikelihood);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + ySlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + zSlope * obsIncrements[i]; }
        double xSlope = regressionSlope(output.yValues, output.xValues);
        zSlope = regressionSlope(output.yValues, output.zValues);
        obsIncrements = this.getObservationIncrements(output.yValues, this.yLikelihood);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + xSlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + zSlope * obsIncrements[i]; }
        xSlope = regressionSlope(output.zValues, output.xValues);
        ySlope = regressionSlope(output.zValues, output.yValues);
        obsIncrements = this.getObservationIncrements(output.zValues, this.zLikelihood);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + xSlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + ySlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        return output;
    }

    /** 
     * Gets observation increments for one variable after finding the posterior
     */
    double[] getObservationIncrements(double[] priorValues, double[] likelihood) {
        //Sort ensemble values and associate them with likelihood
        //I've left a bunch of commented writeln statements here for debugging; will remove later
        assert(likelihood.length == priorValues.length);
        uint[] indexList;
        foreach(i; 0..priorValues.length) {
            indexList ~= cast(uint)i;
        }
        Tuple!(double[], uint[]) sortedLists = indexSort!(double, uint)(priorValues, indexList);
        double[] sortedPrior = sortedLists[0];
        uint[] indices = sortedLists[1];
        //Assign a probability density to each bin
        double[] likelihoodDensity;
        foreach(i; 0..(sortedPrior.length - 1)) {
            likelihoodDensity ~= (likelihood[indices[i + 1]] + likelihood[indices[i]]) / 2;
        }
        //Compute partial Gaussian kernels for tails of prior
        immutable double distanceForUnitSpread = -1 * weightedNormInverse(1, 0, 1, 1.0 / (priorValues.length + 1));
        immutable double leftMean = sortedPrior[0] + distanceForUnitSpread * standardDeviation(sortedPrior, 1);
        immutable double leftStandardDeviation = standardDeviation(sortedPrior, 1);
        immutable double rightMean = sortedPrior[$ - 1] - distanceForUnitSpread * standardDeviation(sortedPrior, 1);
        immutable double rightStandardDeviation = standardDeviation(sortedPrior, 1); 
        //Assume flat tails in likelihood. TODO: Allow for Gaussian tails (should only be relevant for nearly Gaussian likelihoods.)
        immutable double leftProductWeight = likelihood[indices[0]];    
        immutable double rightProductWeight = likelihood[indices[$ - 1]];
        //Assign mass to each bin
        double[] mass = [leftProductWeight / (sortedPrior.length + 1)];
        foreach(i; 0..likelihoodDensity.length) {
            mass ~= likelihoodDensity[i] / (sortedPrior.length + 1);
        }
        mass ~= rightProductWeight / (sortedPrior.length + 1);
        //Get height and normalize mass for trapezoidal, define posterior
        double[] height;
        foreach(i; 0..(sortedPrior.length - 1)) {
            height ~= sortedPrior[i + 1] == sortedPrior[i] ? -1 : 1 / ((sortedPrior.length + 1) * (sortedPrior[i + 1] - sortedPrior[i]));
        }
        double massSum = mass.sum;
        if(massSum == 0) mass[] = 1 / mass.length;
        mass = mass.map!(a => a / massSum).array;
        //Get weight for normalized partial Gaussian tails
        double leftAmp = leftProductWeight / massSum;
        double rightAmp = rightProductWeight / massSum;
        //Get discretely defined CDF of the posterior
        double[] cumulativeMass = [0.0];
        foreach(i; 0..mass.length) {
            cumulativeMass ~= cumulativeMass[i] + mass[i];
        }
        //writeln("LikelihoodDensity: ", likelihoodDensity);
        //writeln("Mass: ", mass);
        //writeln("Cumulative mass: ", cumulativeMass);
        //writeln("priorValues: ", priorValues);
        //writeln("sorted prior: ", sortedPrior);
        //writeln("indices: ", indices);
        //Searching for boxes in which to put each ensemble member
        //TODO: lowest box
        double[] posteriorPoints;
        foreach(i; 0..sortedPrior.length) {
            bool found = false;
            double passedMass = (i + 1.0) / (sortedPrior.length + 1);
            if(passedMass < cumulativeMass[1]) {
                posteriorPoints ~= weightedNormInverse(leftAmp, leftMean, leftStandardDeviation, passedMass);
                assert(!isNaN(cast(float)weightedNormInverse(leftAmp, leftMean, leftStandardDeviation, passedMass)), "RHF(149): Weighted Norm Inverse is NaN");
                found = true;
            } else if(passedMass > cumulativeMass[$ - 2]) {
                posteriorPoints ~= 2 * rightMean - weightedNormInverse(rightAmp, rightMean, rightStandardDeviation, 1 - passedMass);
                assert(rightAmp != 0, "Right alpha is 0");
                assert(!isNaN(cast(float)(2 * rightMean - weightedNormInverse(rightAmp, rightMean, rightStandardDeviation, 1 - passedMass))), "RHF(153): Weighted Norm Inverse is NaN");
                found = true;
            } else {
                foreach(cumulativeMassIndex; 2..cumulativeMass.length) {
                    if(passedMass >= cumulativeMass[cumulativeMassIndex - 1] && passedMass <= cumulativeMass[cumulativeMassIndex] && !found) {
                        if(this.rectangularQuadrature 
                            || (height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 2]])
                            .approxEqual(height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 1]])) {
                            posteriorPoints ~= sortedPrior[cumulativeMassIndex - 2] + (passedMass - cumulativeMass[cumulativeMassIndex - 1]) / (cumulativeMass[cumulativeMassIndex] - cumulativeMass[cumulativeMassIndex - 1]) * (sortedPrior[cumulativeMassIndex - 1] - sortedPrior[cumulativeMassIndex - 2]);
                            found = true;
                        } else {
                            //We're using trapezoidal quadrature to get the new point. 
                            //box is index of cumulative mass, box - 1 is ensemble point
                            double leftHeight = height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 2]] / massSum;
                            double rightHeight = height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 1]] / massSum;
                            //Solve a quadratic to get its roots. Quadratic is the integral of the line that forms the top of a trapezoidal bin
                            double a = 0.5 * (rightHeight - leftHeight) / (sortedPrior[cumulativeMassIndex - 1] - sortedPrior[cumulativeMassIndex - 2]);
                            double b = leftHeight;
                            double c = cumulativeMass[cumulativeMassIndex - 1] - passedMass;
                            Tuple!(double, double) roots = solveQuadratic(a, b, c);
                            double root1 = roots[0];
                            double root2 = roots[1];
                            root1 += sortedPrior[cumulativeMassIndex - 2];
                            root2 += sortedPrior[cumulativeMassIndex - 2];
                            if(root1 >= sortedPrior[cumulativeMassIndex - 2] && root1 <= sortedPrior[cumulativeMassIndex - 1]) {
                                posteriorPoints ~= root1;
                                found = true;
                            } else if(root2 >= sortedPrior[cumulativeMassIndex - 1] && root2 <= sortedPrior[cumulativeMassIndex]) {
                                posteriorPoints ~= root2;
                                found = true;
                            } else {
                                writeln("Lower x value: ", sortedPrior[cumulativeMassIndex - 2], " upper x value: ", sortedPrior[cumulativeMassIndex - 1], " upper x value 2: ", sortedPrior[cumulativeMassIndex]);
                                writeln("left height: ", leftHeight, " right height: ", rightHeight);
                                writeln("a: ", a, " b: ", b, " c: ", c);
                                writeln("root1: ", root1, " root2: ", root2);
                                writeln("ensemble point index: ", cumulativeMassIndex - 2);
                                writeln("Cumulative mass lower point: ", cumulativeMass[cumulativeMassIndex - 1], " cumulative mass higher point: ", cumulativeMass[cumulativeMassIndex]);
                                writeln("Passed mass: ", passedMass);
                                assert(0, "Failed to get proper roots for trapezoidal quadrature.");
                            }
                        }
                    } 
                    /*else {
                        import std.stdio;
                        writeln(passedMass);
                        writeln(cumulativeMass[cumulativeMassIndex - 1], " ", cumulativeMass[cumulativeMassIndex]);
                        assert(0, "");
                    }*/
                }
            }
        }
        //writeln("Prior values: ", priorValues);
        //writeln("Posterior points: ", posteriorPoints);
        double[] observationIncrements;
        assert(posteriorPoints.filter!(a => isNaN(cast(float) a)).array.length == 0, "RHF(205): NaN in posteriorPoints");
        foreach(i; 0..sortedPrior.length) {
            observationIncrements ~= 0;
        }
        foreach(i; 0..sortedPrior.length) {
            assert(indices[i] < observationIncrements.length);
            assert(posteriorPoints.length == sortedPrior.length);
            observationIncrements[indices[i]] = posteriorPoints[i] - sortedPrior[i];
        }
        return cast(double[])observationIncrements;
    }  

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: RHF");
    RHF rhf1 = new RHF(new Likelihood([0, 0, 0.25, 0.5, 1, 0.5, 0.25, 0], [0, 0, 0.25, 0.5, 1, 0.5, 0.25, 0], [0, 0, 0.25, 0.5, 1, 0.5, 0.25, 0]), true);
    Ensemble test1 = new Ensemble(
        [1, 2, 3, 4, 5, 6 ,7, 8], 
        [1, 2, 3, 4, 5, 6 ,7, 8], 
        [1, 2, 3, 4, 5, 6 ,7, 8]
    );
    /*RHF rhf2 = new RHF([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], false);*/
    /*Ensemble test2 = new Ensemble(
        [0, 0.1, 0.2, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 1],
        [0, 0.1, 0.2, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 1],
        [0, 0.1, 0.2, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 1]
    );*/
    Ensemble result1 = rhf1(test1);
    //Ensemble result2 = rhf2(test2);
    writeln("Ensemble with mean ", test1.eMean, 
    " and standard deviation ", test1.eStandardDeviation, 
    " assimilated to <0, 0, 0> with likelihood standard deviation <1, 1, 1> returns an ensemble with mean ",
    result1.eMean, " and standard deviation ", result1.eStandardDeviation);
    /*writeln("Ensemble with mean ", test2.eMean, 
    " and standard deviation ", test2.eStandardDeviation, 
    " assimilated to <0, 0, 0> with likelihood standard deviation <1, 1, 1> returns an ensemble with mean ",
    result2.eMean, " and standard deviation ", result2.eStandardDeviation);*/
    writeln("RHF observation increments on [1, 3, 2] and [0, 1, 2] are: ", rhf1.getObservationIncrements([1, 3, 2], [0, 1, 2]));
    writeln("RHF observation increments on [0, 1, 2] and [1, 1, 1] are: ", rhf1.getObservationIncrements([0, 1, 2], [1, 1, 2]));

}