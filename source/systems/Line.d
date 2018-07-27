module systems.Line;

import std.conv;
import math.Vector;
import systems.System;

/** 
 * An n-dimensional linear system
 */
class Line(uint dim) : System!dim {
    
    Vector!(double, dim) slope;

    this(Vector!(double, dim) slope = new Vector!(double, dim)(1)) {
        this.slope = slope;
    }

    override Vector!(double, dim) opCall(Vector!(double, dim) state) {
        return new Vector!(double, dim)(this.slope);
    }

    override string toString() {
        return "Linear System (slope = " ~ this.slope.toString ~ ")";
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Line");
    Line!3 line = new Line!3();
    writeln("Application of Line to (1, 1, 1) with parameters (1, 1, 1) results in velocity ", line(new Vector!(double, 3)(1, 1, 1)));

}