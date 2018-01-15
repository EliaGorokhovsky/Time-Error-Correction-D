module data.Vector;

import std.conv;
import mir.ndslice;

/**
 * A three-dimensional vector
 * May represent a position vector (point) or a velocity vector (output of a differentially defined system)
 * TODO: n-dimensional vector support
 */
struct Vector {
    
    double x;
    double y;
    double z;

    /**
     * Gets a string representing the vector
     */
    @property string toString() {
        return "<" ~ this.x.to!string ~ "," ~ this.y.to!string ~ "," ~ this.z.to!string ~ ">";
    }

    /**
     * Gets the components of the vector as an array
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
     * Handles vector addition and subtraction
     */
    Vector opBinary(string op)(Vector other) {
        static if (op == "+") return Vector(this.x + other.x, this.y + other.y, this.z + other.z);
        else static if (op == "-") return Vector(this.x - other.x, this.y - other.y, this.z - other.z);
        else return Vector(0, 0, 0);
    }

    /**
     * Handles scalar addition, subtraction, multiplication and division
     */
    Vector opBinary(string op)(double scalar) {
        static if (op == "*") return Vector(this.x * scalar, this.y * scalar, this.z * scalar);
        else static if (op == "/") return Vector(this.x / scalar, this.y / scalar, this.z / scalar);
        else return Vector(0, 0, 0);
    }
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Vector");
    writeln(Vector(1, 2, 3), " + ", Vector(3, 2, 1), " = ", Vector(1, 2, 3) + Vector(3, 2, 1));
    writeln(Vector(1, 2, 3), " x 2 = ", Vector(1, 2, 3) * 2);

}