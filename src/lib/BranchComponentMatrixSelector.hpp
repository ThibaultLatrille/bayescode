#pragma once

#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchSiteSelector.hpp"
#include "SubMatrix.hpp"

/**
 * \brief A BranchSiteSelector
 *
 * The substitution process varies across sites and across branches. In addition,
 * each site is allocated to one of the K component Accordingly, SubMatrixSelector returns, for a
 * given branch and a given site, a pointer to the substitution matrix corresponding to that site,
 * and under the component being active for that site.
 */
template <class T>
class BranchComponentMatrixSelector : public BranchSiteSelector<T> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices and
    //! branch allocation system
    BranchComponentMatrixSelector(const BidimSelector<T> *inmatrixbidimarray, const Selector<int> *inalloc,
        const Tree &intree)
        : matrixbidimarray(inmatrixbidimarray), alloc(inalloc), tree{intree} {}

    int GetSize() const override { return matrixbidimarray->GetNcol(); }

    const T &GetVal(Tree::BranchIndex branch, int site) const override {
        return matrixbidimarray->GetVal(branch, alloc->GetVal(site));
    }

    //! return underlying tree
    const Tree &GetTree() const override { return tree; };

  private:
    const BidimSelector<T> *matrixbidimarray;
    const Selector<int> *alloc;
    const Tree &tree;
};

/**
 * \brief A Selector of substitution matrices at the root
 *
 * The substitution process varies across sites, each site is allocated to one of the K component.
 * Accordingly, SubMatrixSelector returns, for a given site, a pointer to the substitution matrix corresponding to that site,
 * and under the component being active for that site.
 */

template <class T>
class RootComponentMatrixSelector : public Selector<T> {
  public:
    //! constructor, parameterized by bidim array of substitution matrices
    RootComponentMatrixSelector(const Selector<T> *inmatrixbidimarray, const Selector<int> *inalloc)
        : matrixarray(inmatrixbidimarray), alloc{inalloc} {}

    ~RootComponentMatrixSelector() = default;

    int GetSize() const override { return matrixarray->GetSize(); }

    const T &GetVal(int site) const override { return matrixarray->GetVal(alloc->GetVal(site)); }

  private:
    const Selector<T> *matrixarray;
    const Selector<int> *alloc;
};
