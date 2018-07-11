/**
 * Segment
 */
module math.Segment;

import std.math;
import std.parallelism;
import math.AxisAlignedBoundingBox;
import math.Vector;

/**
 * A segment class that is defined by two 
 */
class Segment(T, uint dimensions) {

    Vector!(T, dimensions) initial; ///The initial point (the vector points to the initial point)
    Vector!(T, dimensions) terminal; ///The terminal point (the vector points to the terminal point)

    /**
     * Returns the direction of the segment from initial to terminal
     */
    @property Vector!(T, dimensions) direction() {
        return terminal - initial;
    }

    /**
     * Constructor for a segment;
     * Takes in the initial point and the terminal point
     */
    this(Vector!(T, dimensions) initial, Vector!(T, dimensions) terminal) {
        this.initial = initial;
        this.terminal = terminal;
    }

    /**
     * Returns whether this segment contains the other point
     * Checks that a segment from origin to point has magnitude between initial and terminal's magnitude
     * Also checks that when point magnitude is 1 and direction magnitude is 1, they have the same or they have negated components
     */
    bool contains(T)(Vector!(T, dimensions) point) {
        immutable pointMag = point.magnitude;
        immutable initialMag = this.initial.magnitude;
        immutable terminalMag = this.terminal.magnitude;
        if (!(pointMag > initialMag && pointMag < terminalMag)
                && !(pointMag < initialMag && pointMag > terminalMag)) {
            return false;
        }
        Vector!(T, dimensions) pointCopy = new Vector!(T, dimensions)(point);
        Vector!(T, dimensions) direcCopy = this.direction;
        pointCopy.magnitude = 1;
        direcCopy.magnitude = 1;
        return pointCopy == direcCopy || pointCopy == -direcCopy;
    }

}

/**
 * TODO: Returns where two segments intersect; if they don't intersect, returns null; needs rref of matrix to be completed first
 */
Vector!(T, dimensions) interesction(T, dimensions)(Segment!(T,
        dimensions) first, Segment!(T, dimensions) second) {
    immutable dot = first.direction.dot(second.direction);
    if (!intersects(new AxisAlignedBoundingBox(first.initial.first.terminal),
            new AxisAlignedBoundingBox(first.initial.first.terminal)) || abs(dot).approxEquals(1)) {
        return null;
    }
    return null;
}
