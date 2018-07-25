module experiment.error.UniformTimeError;

import std.algorithm;
import std.array;
import std.math;
import std.range;
import mir.random;
import mir.random.variable;
import data.Timeseries;
import experiment.error.ErrorGenerator;
import integrators.Integrator;

/**
 * Gets error with a uniform time permutation from a single point
 */
class UniformTimeError(uint dim) : ErrorGenerator!dim {

    Vector!(double, dim) error;
    double timeError;
    Integrator integrator;
    Timeseries!(Vector!(double, dim)) truth;
    double firstTime;
    Random* gen;

    this(double timeError, Vector!(double, dim) error, Timeseries!(Vector!(double, dim)) truth, Integrator integrator, Random* gen) {
        this.timeError = timeError;
        this.error = error;
        this.truth = truth;
        this.integrator = integrator;
        this.firstTime = truth.times[0];
        this.gen = gen;
    }

    /**
     * Gets a point observation at a given time by first getting a time from a uniform distribution and then a position
     */
    override Vector!(double, dim) generate(double time) {
        auto newTime = UniformVariable!double(time - this.timeError, time + this.timeError);
        Vector!(double, dim) base = this.truth.value(time, this.integrator);
        double[dim] vars = iota(0, dim, 1).map!(a => NormalVariable!double(base[a], this.error[a])(*this.gen)).array;
        return new Vector!(double, dim)(vars);
    }

}
