module integrators.Integrator;

import data.Vector;
import systems.System;

/**
 * An overarching definition of what makes an integrator
 * Each integrator needs to be able to return 
 */
abstract class Integrator {

    System slope;   ///Should store a function to return slope from position

}