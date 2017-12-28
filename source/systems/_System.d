module systems.System;

import data.Ensemble;
import data.Vector;

/**
 * An overarching definition of what makes a dynamical system
 * Defines what systems need to implement
 */
abstract class System {

    alias slope this;

    Vector slope();     ///The system's defining equations; takes in a position vector and returns a velocity vector for use in an integrator

}