module systems.TestSystem;

import std.conv;
import math.Vector;
import systems.System;

/** 
 * An n-dimensional linear system
 */
class Line(uint dim) : System!dim {
    
    Vector!(double, dim) slope;

    this(Vector!(double, dim) slope = Vector!(double, dim)(1)) {
        this.slope = slope;
    }

    override Vector opCall(Vector state) {
        return this.slope;
    }

    override string toString() {
        return "Test System (slope = " ~ this.slope.toString ~ ")";
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Line");
    Line!3 line = new Line!3();
    writeln("Application of Line to (1, 1, 1) with parameters (1, 1, 1) results in velocity ", line(new Vector!(double, 3)(1, 1, 1)));

}