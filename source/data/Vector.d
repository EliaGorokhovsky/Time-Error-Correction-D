/**
 * My own vector implementation
 * will clash with the d2d vector implementation
 * This is an exclusively 3-dimensional vector with lightweight implementation
 * Handles magnitude and binary operations only
 * Is used to denote position, velocity, and error in the experiment
 */
module data.Vector;

import std.conv; //Used to convert components to strings for string representation
import std.math; //Used for exponents to calculate magnitude

/**
 * A three-dimensional vector
 * May represent a position vector (point) or a velocity vector (output of a differentially defined system)
 * TODO: n-dimensional vector support
 */
struct Vector {
    
    //The three vector components
    double x;
    double y;
    double z;

    /**
     * Gets a string representing the vector
     * returns <x, y, z>
     */
    @property string toString() {
        return "<" ~ this.x.to!string ~ "," ~ this.y.to!string ~ "," ~ this.z.to!string ~ ">";
    }

    /**
     * Gets the components of the vector as an array
     * returns [x, y, z]
     */
    @property double[] handle() {
        return [this.x, this.y, this.z];
    }

    /**
     * Gets the magnitude of the vector or, if used as a point, distance from origin
     */
    @property double magnitude() {
        return(sqrt(this.x.pow(2) + this.y.pow(2) + this.z.pow(2)));
    }

    /**
     * Handles vector operations
     * Vectors are added, subtracted, multiplied, etc. elementwise
     */
    Vector opBinary(string op)(Vector other) {
        mixin("return Vector(this.x " ~ op ~ " other.x, this.y " ~ op ~ " other.y, this.z " ~ op ~ " other.z);");
    }

    /**
     * Handles scalar operations
     * The operation is applied to every element
     */
    Vector opBinary(string op)(double scalar) {
        mixin("return Vector(this.x " ~ op ~ " scalar, this.y " ~ op ~ " scalar, this.z " ~ op ~ " scalar);");
    }
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Vector");
    writeln(Vector(1, 2, 3), " + ", Vector(3, 2, 1), " = ", Vector(1, 2, 3) + Vector(3, 2, 1));
    writeln(Vector(1, 2, 3), " x 2 = ", Vector(1, 2, 3) * 2);

}