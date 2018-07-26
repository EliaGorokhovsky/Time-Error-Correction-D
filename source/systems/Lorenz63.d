module systems.Lorenz63;

import std.conv;
import math.Vector;
import systems.System;

/** 
 * A three-dimensional hydrodynamics toy model defined by Lorenz (1963)
 * Creates a butterfly attractor in phase space
 */
class Lorenz63 : System!3 {

    double rho;
    double sigma;
    double beta;
    
    this(double rho = 28, double sigma = 10, double beta = 8.0/3.0) {
        this.rho = rho;
        this.sigma = sigma;
        this.beta = beta;
    }

    override Vector!(double, 3) opCall(Vector!(double, 3) state) {
        return new Vector!(double, 3)(
            [this.sigma * (state.y - state.x),
            state.x * (this.rho - state.z) - state.y,
            state.x * state.y - this.beta * state.z]
        );
    }

    override string toString() {
        return "L63(rho = " ~ this.rho.to!string ~ ", sigma = " ~ this.sigma.to!string ~ ", beta = " ~ this.beta.to!string ~ ")";
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Lorenz63");
    Lorenz63 L63 = new Lorenz63();
    writeln("Application of L63 to (1, 1, 1) with parameters (1, 1, 1) results in velocity ", L63(new Vector!(double, 3)(1, 1, 1)));

}