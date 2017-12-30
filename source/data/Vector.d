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
}