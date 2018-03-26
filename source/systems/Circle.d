module systems.Circle;

import data.Vector;
import systems.System;
/** 
 * A three-dimensional non-chaotic system that produces a circle
 */
class Circle : System {
    
    double xMod;
    double yMod;

    this(double xMod = 1, double yMod = 1) {
        this.xMod = xMod;
        this.yMod = yMod;
    }

    override Vector opCall(Vector state) {
        return Vector(
            -state.y,
            state.x,
            1
        );
    }

}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Circle");
    Circle circle = new Circle();
    writeln("Application of Circle to (1, 1, 1) with parameters (1, 1, 1) results in velocity ", circle(Vector(1, 1, 1)));

}