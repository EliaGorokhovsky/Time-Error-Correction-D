module experiment.error.GaussianTimeError;

import std.algorithm;
import std.math;
import mir.random;
import mir.random.variable;
import data.Timeseries;
import data.Vector;
import experiment.error.ErrorGenerator;
import integrators.Integrator;

/**
 * Gets error with a Gaussian permutation from a single point
 */
class GaussianTimeError : ErrorGenerator {

    Vector error;
    double timeError;
    Integrator integrator;
    Timeseries!Vector truth;
    double firstTime;
    Random* gen;

    this(double timeError, Vector error, Timeseries!Vector truth, Integrator integrator, Random* gen) {
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
    override Vector generate(double time) {
        auto newTime = NormalVariable!double(time, this.timeError);
        Vector base = this.truth.value(clamp(newTime(*this.gen), 0, 1000000), this.integrator);
        auto normalX = NormalVariable!double(base.x, this.error.x);
        auto normalY = NormalVariable!double(base.y, this.error.y);
        auto normalZ = NormalVariable!double(base.z, this.error.z);
        return Vector(normalX(*this.gen), normalY(*this.gen), normalZ(*this.gen));
    }

}
