module integrators.Integrator;

import data.Ensemble;
import data.Vector;
import systems.System;

/**
 * An overarching definition of what makes an integrator
 * Each integrator needs to be able to return a state from a state and a slope
 */
abstract class Integrator {

    System slope;   ///Should store a function to return slope from position

    Vector integrate(Vector state, double dt);
    Ensemble integrateEnsemble(Ensemble ensemble, double dt);

}