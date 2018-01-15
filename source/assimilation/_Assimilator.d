module assimilation.Assimilator;

import assimilation.likelihood.Likelihood;
import data.Ensemble;
import data.Vector;

/**
 * An overarching definition of what makes a data assimilation method
 * Defines the function of an assimilation step
 */
abstract class Assimilator {

    Ensemble opCall(Ensemble prior);    ///All data assimilation methods should have functionality for each assimilation step
    void setLikelihood(Likelihood likelihood);

}