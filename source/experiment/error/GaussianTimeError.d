module experiment.error.GaussianTimeError;

import std.algorithm;
import std.array;
import std.math;
import std.range;
import std.file;
import std.stdio;
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
class GaussianTimeError(uint dim) : ErrorGenerator!dim {

    Vector!(double, dim) error;
    double timeError;
    Integrator!dim integrator;
    Timeseries!(Vector!(double, dim)) truth;
    double firstTime;
    Random* gen;

    this(double timeError, Vector!(double, dim) error, Timeseries!(Vector!(double, dim)) truth, Integrator!dim integrator, Random* gen) {
        this.timeError = timeError;
        this.error = error;
        this.truth = truth;
        this.integrator = integrator;
        this.firstTime = truth.times[0];
        this.gen = gen;
    }

    /**
     * Gets a point observation at a given time by first getting a time from a normal distribution and then a position
     */
    override Tuple!(double, Vector!(double, dim)) generate(double time) {
        auto newTime = NormalVariable!double(time, this.timeError)(*this.gen);
        Vector!(double, dim) base = this.truth.value(clamp(newTime, 0, 10000000000), this.integrator);
        double[dim] vars = iota(0, dim, 1).map!(a => NormalVariable!double(base[a], this.error[a])(*this.gen)).array;
        return tuple(newTime, new Vector!(double, dim)(vars));
    }

}
