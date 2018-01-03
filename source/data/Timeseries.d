module data.Timeseries;

/**
 * Essentially an array of something
 * Describes the state of an object over a number of steps
 * Contains a variety of useful analysis methods
 */
class Timeseries(T) {

    T[] members; ///The list of the object's states 

    alias members this;

    /**
     * Constructs a timeseries from an existing list of states
     */
    this(T[] members = null) {
        this.members = members;
    }

    /**
     * Adds an element to the timeseries
     */
    void add(T state) {
        this.members ~= state;
    }

}