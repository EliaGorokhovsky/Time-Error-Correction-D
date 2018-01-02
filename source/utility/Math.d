module utility.Math;

import std.math;
import std.typecons;

Tuple!(double, double) solveQuadratic(double a, double b, double c) {
    return tuple((-b + sqrt(b * b - 4 * a * c)) / (2 * a), (-b - sqrt(b * b - 4 * a * c)) / (2 * a));
}