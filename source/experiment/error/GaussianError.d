module experiment.error.GaussianError;

import std.algorithm;
import mir.random;
import mir.random.variable;
import data.Timeseries;
import data.Vector;
import experiment.error.ErrorGenerator;
import integrators.Integrator;

/**
 * Gets error with a Gaussian permutation from a single point
 */
class GaussianError : ErrorGenerator {

    Vector error;
    Integrator integrator;
    Timeseries!Vector truth;

    this(Vector error, Timeseries!Vector truth, Integrator integrator) {
        this.error = error;
        this.truth = truth;
        this.integrator = integrator;
    }

    /**
     * Gets a point observation at a given time
     */
    override Vector generate(double time) {
        Vector base = this.truth.value(time, this.integrator);
        auto gen = Random(unpredictableSeed);
        auto normalX = NormalVariable!double(base.x, this.error.x);
        auto normalY = NormalVariable!double(base.y, this.error.y);
        auto normalZ = NormalVariable!double(base.z, this.error.z);
        return Vector(normalX(gen), normalY(gen), normalZ(gen));
    }

}
