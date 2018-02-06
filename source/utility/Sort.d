module utility.Sort;

import std.algorithm;
import std.array;
import std.math;
import std.typecons;

/**
 * Sort a list by another's values
 */
Tuple!(T[], U[]) indexSort(T, U)(T[] keyList, U[] valueList) {
    assert(keyList.length == valueList.length);
    Tuple!(T, U)[] combinedList;
    foreach(i; 0..keyList.length) {
        combinedList ~= Tuple!(T, U)(keyList[i], valueList[i]);
    }
    assert(combinedList.filter!(a => isNaN(cast(float) a[0])).array.length == 0, "Sort(36): NaN in list");
    combinedList = combinedList.sort!((a, b) => a[0] < b[0])().array;
    T[] sortedKeyList;
    U[] sortedValueList;
    foreach(i; 0..keyList.length) {
        sortedKeyList ~= combinedList[i][0];
        sortedValueList ~= combinedList[i][1];
    }
    return tuple(sortedKeyList, sortedValueList);
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Sort");
    writeln("Sort: [3, 2, 1, 4] [1, 2, 3, 4] => ", indexSort([3, 2, 1, 4], [1, 2, 3, 4]));

}