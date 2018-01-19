module experiment.error.GaussianTimeError;

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

    this(double timeError, Vector error, Timeseries!Vector truth, Integrator integrator) {
        this.error = error;
        this.truth = truth;
        this.integrator = integrator;
    }

    /**
     * Gets a point observation at a given time by firt getting a time from a normal distribution and then a position
     */
    override Vector generate(double time) {
        auto gen = Random(unpredictableSeed);
        auto newTime = NormalVariable!double(time, this.timeError);
        Vector base = this.truth.value(newTime(gen), this.integrator);
        auto normalX = NormalVariable!double(base.x, this.error.x);
        auto normalY = NormalVariable!double(base.y, this.error.y);
        auto normalZ = NormalVariable!double(base.z, this.error.z);
        return Vector(normalX(gen), normalY(gen), normalZ(gen));
    }

}