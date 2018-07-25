module systems.System;

import data.Ensemble;
import math.Vector;

/**
 * An overarching definition of what makes a dynamical system
 * Defines what systems need to implement
 */
abstract class System(uint dim) {

    Vector!(double, dim) opCall(Vector!(double, dim) state);     ///The system's defining equations; takes in a position vector and returns a velocity vector for use in an integrator

}