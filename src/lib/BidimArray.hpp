#pragma once

#include <cassert>
#include "Array.hpp"

/**
 * \brief An interface for an abstract indexed bidim-array of constant
 * references over objects of type T
 *
 * This abstract class is meant as a general read-only interface
 * returning a reference over a T object for any pair of indices (through the
 * GetVal(int i, int j) method). As an example, DiffSelSparseModel implements a
 * bidim array of site- and condition-specific fitness vectors over the 20 amino
 * acids (fitness vector)
 */

template <class T>
class BidimSelector {
  public:
    virtual ~BidimSelector() {}

    //! return number of rows
    virtual int GetNrow() const = 0;
    //! return number of columns
    virtual int GetNcol() const = 0;
    //! const access to element array
    virtual const T &GetVal(int i, int j) const = 0;

    //! write array into generic output stream
    void ToStream(std::ostream &os) const {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) { os << this->GetVal(i, j) << '\t'; }
        }
    }
};

/**
 * \brief An interface for an abstract indexed array of references over objects
 * of type T
 *
 * Deriving from BidimSelector, BidimArray gives const access through GetVal(int
 * i, int j). In addition, it provides a non-const access to through the ()
 * operator. In practice, this class is (should be?) used only for implementing
 * arrays of all-distinct instances of objects of type T. The SimpleBidimArray
 * class template provides a more explicit implementation (by encapsulating a
 * vector<vector<T> >, see SimpleBidimArray). However, in some cases, it can be
 * useful to implement a BidimArray<T> as a vector<vector<T*> > instead, hence
 * the specific declaration of this intermediate class.
 */

template <class T>
class BidimArray : public BidimSelector<T> {
  public:
    virtual ~BidimArray() {}

    //! non-const access to array element
    virtual T &operator()(int i, int j) = 0;

    //! get array from generic input stream
    void FromStream(std::istream &is) {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) { is >> (*this)(i, j); }
        }
    }
};

/**
 * \brief output stream operator for a BidimSelector<T>
 *
 * in practice, calls the ToStream method of BidimSelector<T>
 */

template <class T>
std::ostream &operator<<(std::ostream &os, const BidimSelector<T> &array) {
    array.ToStream(os);
    return os;
}

/**
 * \brief input stream operator for a BidimArray<T>
 *
 * in practice, calls the FromStream method of BidimArray<T>
 */

template <class T>
std::istream &operator>>(std::istream &is, BidimArray<T> &array) {
    array.FromStream(is);
    return is;
}

/**
 * \brief The 'standard' implementation of a BidimArray<T>, simply as a
 * vector<vector<T> >
 *
 */

template <class T>
class SimpleBidimArray : public BidimArray<T> {
  public:
    //! Constructor with array dimensions and initializer value
    SimpleBidimArray(int innrow, int inncol, const T &initval)
        : nrow(innrow), ncol(inncol), array(innrow, std::vector<T>(inncol, initval)) {}
    virtual ~SimpleBidimArray() {}

    int GetNrow() const override { return nrow; }
    int GetNcol() const override { return ncol; }

    T &operator()(int i, int j) override { return array[i][j]; }
    const T &GetVal(int i, int j) const override { return array[i][j]; }

    //! return a const ref to the vector<T> corresponding to row i
    const std::vector<T> &GetSubArray(int i) const { return array[i]; }

    std::vector<std::vector<T>> &GetArray() { return array; }

    template <class Info>
    void declare_interface(Info info) {
        declare(info, "array", array);
    }

  protected:
    int nrow;
    int ncol;
    std::vector<std::vector<T>> array;
};

/**
 * \brief A Selector that distributes the components of a double mixture (double allocation vectors)
 * over an array of items, through 2 vectors of allocations.
 *
 * This Selector<T> maintains a pointer over a Bidimarray of (say K * L) components and
 * 2 arrays of (say N) integer allocation variables. Then, for any index i,
 * GetVal(i) returns a reference to the component to which item i is allocated (using the double
 * allocation) (i.e. components[alloc_row[i], alloc_col[i]])
 */

template <class T>
class DoubleMixtureSelector : public Selector<T> {
  public:
    //! Constructor takes the array of components and the allocation vector
    DoubleMixtureSelector(const BidimSelector<T> *incomponents, const Selector<int> *inalloc_row,
        const Selector<int> *inalloc_col)
        : components(incomponents), alloc_row(inalloc_row), alloc_col(inalloc_col) {
        assert(alloc_col->GetSize() == alloc_row->GetSize());
    }
    ~DoubleMixtureSelector() {}

    //! return size of array of components
    int GetSize() const override { return alloc_col->GetSize(); }

    //! GetVal(i) (for i in 0..alloc->GetSize()-1) returns a reference to the
    //! component to which item i is allocated (i.e. components[alloc_row[i], alloc_col[i]])
    const T &GetVal(int i) const override {
        return components->GetVal(alloc_row->GetVal(i), alloc_col->GetVal(i));
    }

  private:
    const BidimSelector<T> *components;
    const Selector<int> *alloc_row;
    const Selector<int> *alloc_col;
};
