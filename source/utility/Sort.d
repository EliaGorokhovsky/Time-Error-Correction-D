module utility.Sort;

import std.algorithm;
import std.typecons;

/**
 * A struct to hold a key and a value for use in non-uniquely keyed associative arrays
 */
struct keyedValue(T, U) {

    T key;
    U value;

    int opCmp(ref const keyedValue s) const {  
        return cast(int)(this.key - s.key);
    }

    bool opEquals(const keyedValue s) {
        return this.key == s.key;
    }

}

/**
 * Sort a list by another's values
 */
Tuple!(T[], U[]) indexSort(T, U)(T[] keyList, U[] valueList) {
    assert(keyList.length == valueList.length);
    keyedValue!(T, U)[] combinedList;
    foreach(i; 0..keyList.length) {
        combinedList ~= keyedValue!(T, U)(keyList[i], valueList[i]);
    }
    combinedList.sort();
    T[] sortedKeyList;
    U[] sortedValueList;
    foreach(i; 0..keyList.length) {
        sortedKeyList ~= combinedList[i].key;
        sortedValueList ~= combinedList[i].value;
    }
    return tuple(sortedKeyList, sortedValueList);
}

unittest {

    import std.stdio;

    writeln("\nUNITTEST: Sort");
    writeln("Sort: [3, 2, 1, 4] [1, 2, 3, 4] => ", indexSort([3, 2, 1, 4], [1, 2, 3, 4]));

}