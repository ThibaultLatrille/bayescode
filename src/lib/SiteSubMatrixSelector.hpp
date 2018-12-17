#pragma once

#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchSiteSelector.hpp"
#include "SubMatrix.hpp"

/**
 * \brief A BranchSiteSelector for DiffSelModel and DiffSelSparseModel
 *
 * In DiffSelModel and DiffSelSparseModel, the substitution process varies
 * across sites and across the k=0..K-1 conditions. In addition, each branch is
 * allocated to one of the K conditions
 * Accordingly, SubMatrixSelector returns, for a given branch and a given site,
 * a pointer to the substitution matrix corresponding to that site, and under
 * the condition being active for that branch.
 */
template <class T>
class SiteSubMatrixSelector : public BranchSiteSelector<T> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices and
    //! branch allocation system
    SiteSubMatrixSelector(const BidimSelector<T> *inmatrixbidimarray,
        const BranchSelector<int> *inbranchalloc, const Selector<int> *inalloc)
        : matrixbidimarray(inmatrixbidimarray), branchalloc(inbranchalloc), alloc(inalloc) {}

    int GetSize() const override { return matrixbidimarray->GetNcol(); }

    const T &GetVal(int condition, int site) const override {
        return matrixbidimarray->GetVal(branchalloc->GetVal(condition), alloc->GetVal(site));
    }

    //! return underlying tree
    const Tree &GetTree() const override { return branchalloc->GetTree(); };

  private:
    const BidimSelector<T> *matrixbidimarray;
    const BranchSelector<int> *branchalloc;
    const Selector<int> *alloc;
};

/**
 * \brief A Selector of substitution matrices at the root, for DiffSelModel and
 * DiffSelSparseModel
 *
 * In DiffSelModel and DiffSelSparseModel, the root of the tree is assumed to be
 * under condition 0. Accorddingly, RootSubMatrixSelector returns, for a given
 * site, a pointer to the substitution matrix corresponding to that site, and
 * under condition 0.
 */

template <class T>
class RootSiteSubMatrixSelector : public Selector<T> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices
    RootSiteSubMatrixSelector(
        const BidimSelector<T> *inmatrixbidimarray, const Selector<int> *inalloc)
        : matrixbidimarray(inmatrixbidimarray), alloc{inalloc} {}

    ~RootSiteSubMatrixSelector() = default;

    int GetSize() const override { return matrixbidimarray->GetNcol(); }

    const T &GetVal(int site) const override {
        return matrixbidimarray->GetVal(0, alloc->GetVal(site));
    }

  private:
    const BidimSelector<T> *matrixbidimarray;
    const Selector<int> *alloc;
};
