module integrators.RK4;

import std.algorithm;
import std.array;
import data.Ensemble;
import integrators.Integrator;
import math.Vector;
import systems.System;

/**
 * A fourth-order Runge Kutta method for numerical solution of differential equations
 * TODO: adjust for time-sensitive outputs and things with more or less than 3 variables
 */
class RK4(uint dim) : Integrator!dim {

    this(System!dim slope) {
        this.slope = slope;
    }

    /**
     * The integrator's main method
     * Used for numerical solution of differential equations
     * Accurate to O(dt^5)
     * Returns a state given a state and a change in time
     */
    override Vector!(double, dim) integrate(Vector!(double, dim) state, double dt) {
        Vector!(double, dim) k1 = this.slope(state) * dt;
        Vector!(double, dim) k2 = this.slope(state + k1 * 0.5) * dt;
        Vector!(double, dim) k3 = this.slope(state + k2 * 0.5) * dt;
        Vector!(double, dim) k4 = this.slope(state + k3) * dt;
        return state + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6;
    }

    /**
     * Integrates an ensemble 
     */
    override Ensemble!dim integrateEnsemble(Ensemble!dim ensemble, double dt) {
        return new Ensemble!dim(ensemble.members.map!(a => this.integrate(a, dt)).array);
    }

    /**
     * Integrates a given point to a given time
     * Avoid negative timeDifferences, which don't work as expected
     */
    override Vector!(double, dim) integrateTo(Vector!(double, dim) state, double timeDifference, uint steps) {
        Vector!(double, dim) newState = state;
        double dt = timeDifference / steps;
        foreach(i; 0..steps) {
            newState = this.integrate(newState, dt);
        }
        return newState;
    }

    /**
     * Returns the string representation of this integrator
     */
    override string toString() {
        return "RK4";
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: RK4");
    writeln("Differentially defined system: dx/dt = 1, dy/dt = 2, dz/dt = 3");
    class Test : System!3 {
        override Vector!(double, 3) opCall(Vector!(double, 3) state) { return new Vector!(double, 3)([1, 2, 3]); }
    }
    RK4!3 rk4 = new RK4!3(new Test);
    writeln("Base (0, 0, 0) with dt = 1: " ~ rk4(new Vector!(double, 3)(0), 1).toString);
    writeln("Base (0, 0, 0) with dt = 0.5: " ~ rk4(new Vector!(double, 3)(0), 0.5).toString);
    writeln("Base (1, 2, 3) with dt = 1: " ~ rk4(new Vector!(double, 3)([1, 2, 3]), 1).toString);
    writeln("Base Ensemble((3, 2, 1), (1, 2, 3)) with dt = 1: " ~ rk4(new Ensemble!3([new Vector!(double, 3)([1, 2, 3]), new Vector!(double, 3)([3, 2, 1])]), 1).toString);

}