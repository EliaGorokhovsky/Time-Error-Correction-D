module data.Vector;

import std.conv;

/**
 * A three-dimensional vector
 * May represent a position vector (point) or a velocity vector (output of a differentially defined system)
 * TODO: n-dimensional vector support
 */
struct Vector {
    
    double x;
    double y;
    double z;

    @property string toString() {
        return "(" ~ this.x.to!string ~ "," ~ this.y.to!string ~ "," ~ this.z.to!string ~ ")";
    }

}