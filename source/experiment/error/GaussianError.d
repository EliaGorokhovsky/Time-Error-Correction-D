module experiment.error.GaussianError;

import std.algorithm;
import std.array;
import std.range;
import std.typecons;
import mir.random;
import mir.random.variable;
import data.Timeseries;
import experiment.error.ErrorGenerator;
import integrators.Integrator;
import math.Vector;

/**
 * Gets error with a Gaussian permutation from a single point
 */
class GaussianError(uint dim) : ErrorGenerator!dim {

    Vector!(double, dim) error;
    Integrator!dim integrator;
    Timeseries!(Vector!(double, dim)) truth;
    Random* gen;

    this(Vector!(double, dim) error, Timeseries!(Vector!(double, dim)) truth, Integrator!dim integrator, Random* gen) {
        this.error = error;
        this.truth = truth;
        this.integrator = integrator;
        this.gen = gen;
    }

    /**
     * Gets a point observation at a given time
     */
    override Tuple!(double, Vector!(double, dim)) generate(double time) {
        Vector!(double, dim) base = this.truth.value(time, this.integrator);
        double[dim] vars = iota(0, dim, 1).map!(a => NormalVariable!double(base[a], this.error[a])(*this.gen)).array;
        return tuple(time, new Vector!(double, dim)(vars));
    }

}
