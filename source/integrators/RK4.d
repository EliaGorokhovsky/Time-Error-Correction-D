module integrators.RK4;

import std.algorithm;
import std.array;
import data.Ensemble;
import integrators.Integrator;
import systems.System;

/**
 * A fourth-order Runge Kutta method for numerical solution of differential equations
 * TODO: adjust for time-sensitive outputs and things with more or less than 3 variables
 */
class RK4 : Integrator {

    this(System slope) {
        this.slope = slope;
    }

    /**
     * The integrator's main method
     * Used for numerical solution of differential equations
     * Accurate to O(dt^5)
     * Returns a state given a state and a change in time
     */
    override Vector integrate(Vector state, double dt) {
        double[] k1 = this.slope(state).handle.map!(a => a * dt).array;
        double[] k2 = this.slope(
            Vector(state.x + 0.5 * k1[0], state.y + 0.5 * k1[1], state.z + 0.5 * k1[2])
            ).handle.map!(a => a * dt).array;
        double[] k3 = this.slope(
            Vector(state.x + 0.5 * k2[0], state.y + 0.5 * k2[1], state.z + 0.5 * k2[2])
            ).handle.map!(a => a * dt).array;
        double[] k4 = this.slope(
            Vector(state.x + k3[0], state.y + k3[1], state.z + k3[2])
            ).handle.map!(a => a * dt).array;
        return Vector(
            state.x + k1[0] / 6.0 + k2[0] / 3.0 + k3[0] / 3.0 + k4[0] / 6.0,
            state.y + k1[1] / 6.0 + k2[1] / 3.0 + k3[1] / 3.0 + k4[1] / 6.0,
            state.z + k1[2] / 6.0 + k2[2] / 3.0 + k3[2] / 3.0 + k4[2] / 6.0
        );
    }

    /**
     * Integrates an ensemble 
     */
    override Ensemble integrateEnsemble(Ensemble ensemble, double dt) {
        return new Ensemble(ensemble.members.map!(a => this.integrate(a, dt)).array);
    }

    /**
     * Integrates a given point to a given time
     * Avoid negative timeDifferences, which don't work as expected
     */
    override Vector integrateTo(Vector state, double timeDifference, uint steps) {
        Vector newState = state;
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
    class Test : System {
        override Vector opCall(Vector state) { return Vector(1, 2, 3); }
    }
    RK4 rk4 = new RK4(new Test);
    writeln("Base (0, 0, 0) with dt = 1: " ~ rk4(Vector(0, 0, 0), 1).toString);
    writeln("Base (0, 0, 0) with dt = 0.5: " ~ rk4(Vector(0, 0, 0), 0.5).toString);
    writeln("Base (1, 2, 3) with dt = 1: " ~ rk4(Vector(1, 2, 3), 1).toString);
    writeln("Base Ensemble((3, 2, 1), (1, 2, 3)) with dt = 1: " ~ rk4(new Ensemble([Vector(1, 2, 3), Vector(3, 2, 1)]), 1).toString);

}