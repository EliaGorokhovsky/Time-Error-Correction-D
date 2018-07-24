/**
 * Contains the Rank Histogram Filter
 * Which was written by Dr. Jeff Anderson of NCAR
 * The RHF creates a prior distribution by splitting equal mass in between ensemble members in 1 dimension
 * Then it finds a posterior by pairwise multiplication between likelihood and prior
 * And then spreads the senbmel points so that there is equal mass between them in the posterior
 * A rank histogram is a way to approximate the distribution of a set of points by
 * placing equal mass in between each pair of points to create a histogram
 */
module assimilation.RHF;

import std.algorithm; //Used for map statements
import std.array; //Used to convert map statement outputs to arrays
import std.conv; //Used for type conversion 
import std.mathspecial; //Used for approximate equality between doubles and checking for NaN
import std.range; //Used for iotas (ranges of numbers with steps)
import std.typecons; //Used for tuples for sorting pairs of items
import assimilation.Assimilator; //Parent class
import assimilation.likelihood.Likelihood; //Used for input likelihood
import data.Ensemble; //Used for input and output ensemble
import experiment.Analytics; //Used to check if an array has NaN values
import utility.ArrayStats; //Used for standard deviation calculations for Gaussian tails
import utility.Math; //Used for approximate equality between doubles
import utility.Normal; //Used to get the value for the normal probability density function when computing Gaussian tails
import utility.Regression; //Used to regress observation increments onto state variables
import utility.Sort; //Used to sort associated pairs of items

/** 
 * The Rank Histogram Filter (Anderson 2010) to apply ensemble filtering with discretely defined likelihood distributions
 */
class RHF : Assimilator {

    double[] xLikelihood; ///The probabilities for the x-values of the ensemble
    double[] yLikelihood; ///The probabilities for the y-values of the ensemble
    double[] zLikelihood; ///The probabilities for the z-values of the ensemble
    bool rectangularQuadrature; ///Whether rectangular or trapezoidal quadrature will be used to calculate probabilities with likelihood

    /**
     * Constructs the RHF
     * By default, rectangular quadrature is true
     */
    this(Likelihood likelihood, bool rectangularQuadrature = true) {
        this.xLikelihood = likelihood.xLikelihood;
        this.yLikelihood = likelihood.yLikelihood;
        this.zLikelihood = likelihood.zLikelihood;
        this.rectangularQuadrature = rectangularQuadrature;
    }

    /**
     * The empty constructor for an assimilator
     * We can construct an assimilator without any input, then set the likelihood later
     * if we intend to reuse it
     */
    this(bool rectangularQuadrature = true) {
        this.rectangularQuadrature = rectangularQuadrature;
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
        //Ensure none of the likelihoods have NaN (not a number) in them
        assert(!checkNaN(this.xLikelihood), "NaN in X likelihood");
        assert(!checkNaN(this.xLikelihood), "NaN in Y likelihood");
        assert(!checkNaN(this.xLikelihood), "NaN in Z likelihood");
        assert(!checkNaN(prior.xValues), "NaN in X prior");
        assert(!checkNaN(prior.yValues), "NaN in Y prior");
        assert(!checkNaN(prior.zValues), "NaN in Z prior");
        //Copy the ensemble so as to not change it while calculating
        Ensemble output = prior.copy();
        //Regress x observation increments onto y and z
        double ySlope = regressionSlope(output.xValues, output.yValues);
        double zSlope = regressionSlope(output.xValues, output.zValues);
        double[] obsIncrements = this.getObservationIncrements(output.xValues, this.xLikelihood);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + ySlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + zSlope * obsIncrements[i]; }
        //Regress y observation increments onto x and z
        double xSlope = regressionSlope(output.yValues, output.xValues);
        zSlope = regressionSlope(output.yValues, output.zValues);
        obsIncrements = this.getObservationIncrements(output.yValues, this.yLikelihood);
        foreach(i; 0..obsIncrements.length) { output.members[i].x = output.members[i].x + xSlope * obsIncrements[i]; }
        foreach(i; 0..obsIncrements.length) { output.members[i].y = output.members[i].y + obsIncrements[i]; } //Regression of a variable onto itself returns 1
        foreach(i; 0..obsIncrements.length) { output.members[i].z = output.members[i].z + zSlope * obsIncrements[i]; }
        //Regress z observation increments onto x and y
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
        //Ensure likelihood fits prior
        assert(likelihood.length == priorValues.length);
        //Sort ensemble values and associate them with likelihood
        uint[] indexList = iota(0, priorValues.length.to!uint, 1).array;
        Tuple!(double[], uint[]) sortedLists = indexSort!(double, uint)(priorValues, indexList);
        double[] sortedPrior = sortedLists[0];
        uint[] indices = sortedLists[1];
        //Assign a probability density to each bin
        //This is the area of the likelihood that falls in the bin, approximated as a trapezoid,
        //compared to the prior
        double[] likelihoodDensity;
        foreach(i; 0..(sortedPrior.length - 1)) {
            likelihoodDensity ~= (likelihood[indices[i + 1]] + likelihood[indices[i]]) / 2;
        }
        //Compute partial Gaussian kernels for tails of prior
        //Our prior is defined discretely between the points, but the ends are defined by partial Gaussian kernels
        //Find the distance between the mean and the ends of the prior so that the standard deviation of the tails is 1
        immutable double distanceForUnitSpread = -1 * weightedNormInverse(1, 0, 1, 1.0 / (priorValues.length + 1));
        //Find the mean of the partial Gaussian kernel on the left tail of the distribution
        immutable double leftMean = sortedPrior[0] + distanceForUnitSpread * standardDeviation!1(sortedPrior);
        //Find the standard deviation of the partial Gaussian kernel on the left such that its area is equal to those of the bins
        immutable double leftStandardDeviation = standardDeviation!1(sortedPrior);
        //Find the mean of the right Gaussian kernel
        immutable double rightMean = sortedPrior[$ - 1] - distanceForUnitSpread * standardDeviation!1(sortedPrior);
        //Find the correct standard deviation for the right Gaussian kernels
        immutable double rightStandardDeviation = standardDeviation!1(sortedPrior); 
        //Assume flat tails in likelihood. TODO: Allow for Gaussian tails (should only be relevant for nearly Gaussian likelihoods.)
        //Find the tails of the likelihood distribution
        immutable double leftProductWeight = likelihood[indices[0]];    
        immutable double rightProductWeight = likelihood[indices[$ - 1]];
        //Assign mass to each bin so that the masses in each bin are equal to each other and the tails
        //Then multiply the masses by the likelihood density
        //We can do all of this in one step
        double[] mass = [leftProductWeight / (sortedPrior.length + 1)];
        foreach(i; 0..likelihoodDensity.length) {
            mass ~= likelihoodDensity[i] / (sortedPrior.length + 1);
        }
        mass ~= rightProductWeight / (sortedPrior.length + 1);
        //Get height and normalize mass for trapezoidal, define posterior
        //The height represents the height of the rectangle in each bin
        double[] height;
        foreach(i; 0..(sortedPrior.length - 1)) {
            //If a bin's width is 0 then the height is marked as -1
            height ~= sortedPrior[i + 1] == sortedPrior[i] ? -1 : 1 / ((sortedPrior.length + 1) * (sortedPrior[i + 1] - sortedPrior[i]));
        }
        //Normalize the rank histogram's mass
        double massSum = mass.sum;
        if(massSum == 0) mass[] = 1 / mass.length;
        mass = mass.map!(a => a / massSum).array;
        //Get weight for normalized partial Gaussian tails
        double leftAmp = leftProductWeight / massSum;
        double rightAmp = rightProductWeight / massSum;
        //Get discretely defined CDF (cumulative distribution function) of the posterior
        double[] cumulativeMass = [0.0];
        foreach(i; 0..mass.length) {
            cumulativeMass ~= cumulativeMass[i] + mass[i];
        }
        //Searching for boxes in which to put each ensemble member
        //So that the ensemble members are evenly distributed by mass in the posterior
        //i.e. the mass between each pair of ensemble members is 1 / (n + 1)
        //We can speed this up by keeping track of the lowest checked box so that we don't redo the first few boxes every time
        //but this is not currently implemented
        double[] posteriorPoints;
        foreach(i; 0..sortedPrior.length) {
            bool found = false;
            //Remember how much mass should be before the new position of this ensemble member
            double passedMass = (i + 1.0) / (sortedPrior.length + 1);
            //If the ensemble member should fall within the left tail of the distribution
            //We can use the inverse normal function to find where it goes
            if(passedMass < cumulativeMass[1]) {
                posteriorPoints ~= weightedNormInverse(leftAmp, leftMean, leftStandardDeviation, passedMass);
                assert(!isNaN(cast(float)weightedNormInverse(leftAmp, leftMean, leftStandardDeviation, passedMass)), "RHF(149): Weighted Norm Inverse is NaN");
                found = true;
            //If the ensemble member should fall within the right tail of the distribution
            //We can use the inverse normal function to find where it goes
            } else if(passedMass > cumulativeMass[$ - 2]) {
                //Reflect the distribution about the mean in order to adjust for the tail being the right tail
                posteriorPoints ~= 2 * rightMean - weightedNormInverse(rightAmp, rightMean, rightStandardDeviation, 1 - passedMass);
                assert(rightAmp != 0, "Right alpha is 0");
                assert(!isNaN(cast(float)(2 * rightMean - weightedNormInverse(rightAmp, rightMean, rightStandardDeviation, 1 - passedMass))), "RHF(153): Weighted Norm Inverse is NaN");
                found = true;
            //Otherwise we need to check all the bins
            } else {
                //Check each bin
                foreach(cumulativeMassIndex; 2..cumulativeMass.length) {
                    //If the ensemble point should fall within the bin, then its passedMass should be below the cumulative mass of the appropriate index
                    if(passedMass >= cumulativeMass[cumulativeMassIndex - 1] && passedMass <= cumulativeMass[cumulativeMassIndex] && !found) {
                        //We can approximate this with either rectangular or trapezoidal area approximation
                        //The kind of quadrature relates to how we see the bins
                        if(this.rectangularQuadrature 
                            || (height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 2]])
                            .approxEqual(height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 1]])) {
                            //If we are using rectangles, OR if we are not but the bin is rectangular anyway
                            //Then the area increases linearly so we can just find its slope
                            posteriorPoints ~= sortedPrior[cumulativeMassIndex - 2] + (passedMass - cumulativeMass[cumulativeMassIndex - 1]) / (cumulativeMass[cumulativeMassIndex] - cumulativeMass[cumulativeMassIndex - 1]) * (sortedPrior[cumulativeMassIndex - 1] - sortedPrior[cumulativeMassIndex - 2]);
                            found = true;
                        } else {
                            //Check if the ensemble falls in between bins
                            if(passedMass.approxEqual(cumulativeMass[cumulativeMassIndex - 1])) {
                                posteriorPoints ~= sortedPrior[cumulativeMassIndex - 2];
                                found = true;
                            } else if(passedMass.approxEqual(cumulativeMass[cumulativeMassIndex])) {
                                posteriorPoints ~= sortedPrior[cumulativeMassIndex - 1];
                                found = true;
                            } else {
                                //We're using trapezoidal quadrature to get the new point. 
                                //box is index of cumulative mass, box - 1 is ensemble point
                                //These are the heights of the trapezoid that the mass in the bin looks like
                                double leftHeight = height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 2]] / massSum;
                                double rightHeight = height[cumulativeMassIndex - 2] * likelihood[indices[cumulativeMassIndex - 1]] / massSum;
                                //Solve a quadratic to get its roots. Quadratic is the integral of the line that forms the top of a trapezoidal bin
                                double a = 0.5 * (rightHeight - leftHeight) / (sortedPrior[cumulativeMassIndex - 1] - sortedPrior[cumulativeMassIndex - 2]);
                                double b = leftHeight;
                                double c = cumulativeMass[cumulativeMassIndex - 1] - passedMass;
                                Tuple!(double, double) roots = solveQuadratic(a, b, c);
                                double root1 = roots[0];
                                double root2 = roots[1];
                                //Shift the roots over to the ensemble member
                                root1 += sortedPrior[cumulativeMassIndex - 2];
                                root2 += sortedPrior[cumulativeMassIndex - 2];
                                //Find the root that is within the bounds of the bin
                                if(root1 >= sortedPrior[cumulativeMassIndex - 2] && root1 <= sortedPrior[cumulativeMassIndex - 1] ||
                                root1.approxEqual(sortedPrior[cumulativeMassIndex - 2]) || root1.approxEqual(sortedPrior[cumulativeMassIndex - 1])) {
                                    posteriorPoints ~= root1;
                                    found = true;
                                } else if(root2 >= sortedPrior[cumulativeMassIndex - 1] && root2 <= sortedPrior[cumulativeMassIndex] || 
                                root1.approxEqual(sortedPrior[cumulativeMassIndex - 1]) || root1.approxEqual(sortedPrior[cumulativeMassIndex])) {
                                    posteriorPoints ~= root2;
                                    found = true;
                                } else {
                                    //Only use the below code when debugging; these statements should never be reached
                                    /*
                                    import std.stdio;
                                    writeln("Prior: ", priorValues);
                                    writeln("Lower x value: ", sortedPrior[cumulativeMassIndex - 2], " upper x value: ", sortedPrior[cumulativeMassIndex - 1], " upper x value 2: ", sortedPrior[cumulativeMassIndex]);
                                    writeln("left height: ", leftHeight, " right height: ", rightHeight);
                                    writeln("a: ", a, " b: ", b, " c: ", c);
                                    writeln("root1: ", root1, " root2: ", root2);
                                    writeln("ensemble point index: ", cumulativeMassIndex - 2);
                                    writeln("Cumulative mass lower point: ", cumulativeMass[cumulativeMassIndex - 1], " cumulative mass higher point: ", cumulativeMass[cumulativeMassIndex]);
                                    writeln("Passed mass: ", passedMass);
                                    */
                                    assert(0, "Failed to get proper roots for trapezoidal quadrature.");
                                }
                            }
                        }
                    } 
                }
            }
        }
        //Comput observation increments by finding the difference between the new and old ensembles
        double[] observationIncrements;
        //Ensure that the posterior calculation was sucessful
        assert(!checkNaN(posteriorPoints), "RHF: NaN in posteriorPoints");
        foreach(i; 0..sortedPrior.length) {
            observationIncrements ~= 0;
        }
        //Ensure that the prior ensemble is the same size as the posterior
        assert(posteriorPoints.length == sortedPrior.length, "There are " ~ posteriorPoints.length.to!string ~ " posterior points and " ~ sortedPrior.length.to!string ~ " prior points.");
        //Compute differences pointwise
        foreach(i; 0..sortedPrior.length) {
            assert(indices[i] < observationIncrements.length);
            observationIncrements[indices[i]] = posteriorPoints[i] - sortedPrior[i];
        }
        return cast(double[])observationIncrements;
    }  

    /**
     * What to return when attempting to print this assimilator
     */
    override string toString() {
        return "RHF";
    }

}

unittest {

    import std.stdio;

    //Ensure results are similar to EAKF
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