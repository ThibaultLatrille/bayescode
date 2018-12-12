
#ifndef BRANCHSITEARRAY_H
#define BRANCHSITEARRAY_H

#include "Array.hpp"
#include "BranchArray.hpp"

/**
 * \brief An interface for an abstract bidim array over branches and over sites,
 * of constant references over objects of type T
 *
 * This abstract class is meant as a general read-only interface
 * returning a reference over a T object for any branch/site pair (through the
 * GetVal(int branch, int site) method). Branches are indexed from 0 to
 * GetNbranch()-1, and sites from 0 to GetSize()-1. This class, and its various
 * implementations/specializations, follows the same logic as the classes
 * deriving from Array.
 *
 * A good example is the const BranchSiteSelector<SubMatrix>& member of
 * PhyloProcess. All calculations performed by PhyloProcess do not depend on the
 * specific modulations of the substitution process across sites and branches
 * implied by each model. As a consequence, the only thing required by
 * PhyloProcess is a generic interface specifying which substitution matrix
 * should be used for a given site and over a given branch.
 *
 * In practice, this interface will hide a great variety of behaviors.
 * Some models will assume the same substitution matrix for all sites and over
 * all branches (in which case the interface will be specialized into a
 * BranchHomogeneousSiteHomogeneousSelector), an array of matrices across sites,
 * but homogeneous across branches (BranchHomogeneousSiteHeterogeneousSelector),
 * or conversely, distinct matrices across branches, but identical for all sites
 * (BranchHeterogeneousSiteHomogeneousSelector).
 */

template <class T>
class BranchSiteSelector {
  public:
    virtual ~BranchSiteSelector() {}

    //! return underlying tree
    virtual const Tree &GetTree() const = 0;
    //! return number of branches of the underlying tree
    int GetNbranch() const { return GetTree().nb_branches(); }

    //! return number of sites
    virtual int GetSize() const = 0;
    //! const access for a given branch/site pair
    virtual const T &GetVal(int branch, int site) const = 0;
};

/**
 * \brief A specialization of BranchSiteSelector, returning the same const
 * reference to an object of type T for all sites and all branches
 */

template <class T>
class BranchHomogeneousSiteHomogeneousSelector : public BranchSiteSelector<T> {
  public:
    BranchHomogeneousSiteHomogeneousSelector(const Tree &intree, int insize, const T &invalue)
        : tree(intree), size(insize), value(invalue) {}
    ~BranchHomogeneousSiteHomogeneousSelector() {}

    const Tree &GetTree() const /*override*/ { return tree; }
    int GetSize() const /*override*/ { return size; }
    const T &GetVal(int branch, int site) const /*override*/ { return value; }

  private:
    const Tree &tree;
    int size;
    const T &value;
};

/**
 * \brief A specialization of BranchSiteSelector, returning the same const
 * reference to an object of type T for all branches, for a given site
 */

template <class T>
class BranchHomogeneousSiteHeterogeneousSelector : public BranchSiteSelector<T> {
  public:
    //! constructor is parameterized by a tree and an Array<T> across sites
    BranchHomogeneousSiteHeterogeneousSelector(const Tree &intree, const Selector<T> &inarray)
        : tree(intree), array(inarray) {}
    ~BranchHomogeneousSiteHeterogeneousSelector() {}

    const Tree &GetTree() const /*override*/ { return tree; }
    int GetSize() const /*override*/ { return array.GetSize(); }
    //! effectively returns array[site], for any branch
    const T &GetVal(int branch, int site) const /*override*/ { return array.GetVal(site); }

  private:
    const Tree &tree;
    const Selector<T> &array;
};

/**
 * \brief A specialization of BranchSiteSelector, returning the same const
 * reference to an object of type T for all sites, for a given branch
 */

template <class T>
class BranchHeterogeneousSiteHomogeneousSelector : public BranchSiteSelector<T> {
  public:
    //! constructor is parameterized by a BranchArray<T> and a number of site
    //! (insize)
    BranchHeterogeneousSiteHomogeneousSelector(const BranchSelector<T> &inbrancharray, int insize)
        : brancharray(inbrancharray), size(insize) {}
    ~BranchHeterogeneousSiteHomogeneousSelector() {}

    const Tree &GetTree() const /*override*/ { return brancharray.GetTree(); }
    int GetSize() const /*override*/ { return size; }
    //! effectively returns brancharray[branch], for any site
    const T &GetVal(int branch, int site) const /*override*/ { return brancharray.GetVal(branch); }

  private:
    const BranchSelector<T> &brancharray;
    int size;
};

#endif
