module experiment.error.ErrorGenerator;

import data.Timeseries;
import math.Vector;

/**
 * A class that deals with generation of noise from a particular point
 */
abstract class ErrorGenerator(uint dim) {

    Timeseries!(Vector!(double, dim)) truth;

    Vector!(double, dim) opCall(double time) {
        return this.generate(time);
    }

    Vector!(double, dim) generate(double time);

}