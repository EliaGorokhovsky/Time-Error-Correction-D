module experiment.error.GaussianError;

import mir.random;
import mir.random.variable;
import data.Timeseries;
import data.Vector;
import experiment.error.ErrorGenerator;

/**
 * Gets error with a Gaussian permutation from a single point
 */
class GaussianError : ErrorGenerator {

    Vector error;

    this(Vector error, Timeseries!Vector truth) {
        this.error = error;
    }

    /**
     * Gets a point observation at a given time
     */
    override Vector opCall(double time) {
        Vector base = this.truth.value(time);
        auto gen = Random(unpredictableSeed);
        auto normalX = NormalVariable!double(base.x, this.error.x);
        auto normalY = NormalVariable!double(base.y, this.error.y);
        auto normalZ = NormalVariable!double(base.z, this.error.z);
        return Vector(normalX(gen), normalY(gen), normalZ(gen));
    }

}