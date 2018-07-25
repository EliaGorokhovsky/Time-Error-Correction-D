module integrators.Integrator;

import data.Ensemble;
import math.Vector;
import systems.System;

/**
 * An overarching definition of what makes an integrator
 * Each integrator needs to be able to return a state from a state and a slope
 */
abstract class Integrator(uint dim) {

    System!dim slope;   ///Should store a function to return slope from position

    Vector!(double, dim) opCall(Vector!(double, dim) state, double dt) {
        return this.integrate(state, dt);
    }

    Ensemble!dim opCall(Ensemble!dim ensemble, double dt) {
        return this.integrateEnsemble(ensemble, dt);
    }

    Vector!(double, dim) integrate(Vector!(double, dim) state, double dt);
    Ensemble!dim integrateEnsemble(Ensemble!dim ensemble, double dt);
    Vector!(double, dim) integrateTo(Vector!(double, dim) state, double timeDifference, uint steps);

}