/**
 * This module contains the parent class of all assimilators
 * It defines what an assimilator needs to have to be acceptable
 * An assimilator is a class that takes a likelihood and an ensemble, and
 * returns an ensemble that is updated with the new information
 */
module assimilation.Assimilator;

import assimilation.likelihood.Likelihood; //Used for inputting likelihoods
import data.Ensemble; //Used for input and output ensembles
import data.Vector; //Used for calculations in 3d

/**
 * An overarching definition of what makes a data assimilation method
 * Defines the function of an assimilation step
 */
abstract class Assimilator {

    Ensemble opCall(Ensemble prior);    ///All data assimilation methods should have functionality for each assimilation step
    void setLikelihood(Likelihood likelihood); ///All data assimilation methods should accept likelihoods for assimilation

}