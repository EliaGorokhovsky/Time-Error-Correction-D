module experiment.error.GaussianError;

import std.algorithm;
import mir.random;
import mir.random.variable;
import data.Timeseries;
import experiment.error.ErrorGenerator;
import integrators.Integrator;

/**
 * Gets error with a Gaussian permutation from a single point
 */
class GaussianError : ErrorGenerator {

    Vector error;
    Integrator integrator;
    Timeseries!Vector truth;
    Random* gen;

    this(Vector error, Timeseries!Vector truth, Integrator integrator, Random* gen) {
        this.error = error;
        this.truth = truth;
        this.integrator = integrator;
        this.gen = gen;
    }

    /**
     * Gets a point observation at a given time
     */
    override Vector generate(double time) {
        Vector base = this.truth.value(time, this.integrator);
        auto normalX = NormalVariable!double(base.x, this.error.x);
        auto normalY = NormalVariable!double(base.y, this.error.y);
        auto normalZ = NormalVariable!double(base.z, this.error.z);
        return Vector(normalX(*this.gen), normalY(*this.gen), normalZ(*this.gen));
    }

}
