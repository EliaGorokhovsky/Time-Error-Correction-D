module utility.NDVector;

import std.algorithm;
import std.array;

/**
 * An n-dimensional vector
 * Essentially a dynamic array with vector-like operations
 */
 class NDVector(T) {

    T[] members; ///The components of the vector
    alias members this;

    /**
     * Constructs a new NDVector from an array
     */
    this(T[] array) {
        this.members = array;
    }

    /**
     * NDVector operations are performed element-wise
     */
    NDVector!T opBinary(string op)(T scalar) {
        static if(op == "~") { return new NDVector!T(this.members ~ scalar); }
        else { mixin("return new NDVector!T(this.members.map!(a => a " ~ op ~ " scalar).array);"); }
    }
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: NDVector");
    NDVector!int vector = new NDVector!int([1, 2, 3]);
    writeln("vector = [1, 2, 3]");
    writeln("vector + 3 = ", vector + 3);
    writeln("vector * 3 = ", vector * 3);
    writeln("vector - 3 = ", vector - 3);
    writeln("vector / 3 = ", vector / 3);
    writeln("vector ~ 4 = ", vector ~ 4);

}