module systems.Lorenz63;

import data.Vector;
import systems.System;
/** 
 * A three-dimensional hydrodynamics toy model defined by Lorenz (1963)
 * Creates a butterfly attractor in phase space
 */
class Lorenz63 : System {

    double rho;
    double sigma;
    double beta;
    
    this(double rho = 28, double sigma = 10, double beta = 8.0/3.0) {
        this.rho = rho;
        this.sigma = sigma;
        this.beta = beta;
    }

    override Vector opCall(Vector state) {
        return Vector(
            this.sigma * (state.y - state.x),
            state.x * (this.rho - state.z) - state.y,
            state.x * state.y - this.beta * state.z
        );
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Lorenz63");
    Lorenz63 L63 = new Lorenz63();
    writeln("Application of L63 to (1, 1, 1) with parameters (1, 1, 1) results in velocity ", L63(Vector(1, 1, 1)));

}