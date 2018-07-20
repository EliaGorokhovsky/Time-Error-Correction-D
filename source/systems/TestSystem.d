module systems.TestSystem;

import std.conv;
import data.Vector;
import systems.System;
/** 
 * A three-dimensional non-chaotic system that produces a circle
 */
class TestSystem : System {
    
    double xSlope;
    double ySlope;
    double zSlope;

    this(double xSlope = 5, double ySlope = 10, double zSlope = 15) {
        this.xSlope = xSlope;
        this.ySlope = ySlope;
        this.zSlope = zSlope;
    }

    override Vector opCall(Vector state) {
        return Vector(
            xSlope,
            ySlope,
            zSlope
        );
    }

    override string toString() {
        return "Test System (xSlope = " ~ this.xSlope.to!string ~ ", ySlope = " ~ this.ySlope.to!string ~ ", zSlope = " ~ this.zSlope.to!string ~ ")";
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Circle");
    Circle circle = new Circle();
    writeln("Application of Circle to (1, 1, 1) with parameters (1, 1, 1) results in velocity ", circle(Vector(1, 1, 1)));

}