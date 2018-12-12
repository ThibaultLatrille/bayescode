#pragma once

#include "Array.hpp"
#include "BranchArray.hpp"

/**
 * \brief The product of a BranchSelector<double> and a double
 *
 * takes two arguments: a branch array l_j and a double multiplier a, and
 * returns, for branch j, the product a * l_j.
 *
 */

class BranchProduct : public BranchSelector<double> {
  public:
    //! \brief Constructor with branch array and multiplier
    BranchProduct(const BranchSelector<double> &inbranchval, double inmulval)
        : branchval(inbranchval), mulval(inmulval), array(inbranchval.GetNbranch()) {}
    ~BranchProduct() {}

    virtual const Tree &GetTree() const override { return branchval.GetTree(); }

    virtual const double &GetVal(int index) const override {
        array[index] = mulval * branchval.GetVal(index);
        return array[index];
    }

    //! \brief specifies a new value for the multiplier
    void SetMulVal(double inmulval) { mulval = inmulval; }

  private:
    const BranchSelector<double> &branchval;
    double mulval;
    mutable std::vector<double> array;
};

/**
 * \brief The product of a BranchSelector<double> and a Selector<double>
 *
 * Takes two arguments: a BranchArray l_j and an Array r_i, and
 * returns, for branch j, and item i, the product r_i * l_j.
 *
 * Used in BranchOmegaModel: the mean of omega_gj, for gene g and branch j, is
 * equal to w_g * branchv_j
 */

class BranchProductArray : public Array<BranchProduct> {
  public:
    //! Constructor: parameterized by the branch array and array
    BranchProductArray(const BranchSelector<double> &inbranchval, const Selector<double> &insiteval)
        : branchval(inbranchval),
          siteval(insiteval),
          array(insiteval.GetSize(), (BranchProduct *)0) {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene] = new BranchProduct(branchval, siteval.GetVal(gene));
        }
    }

    ~BranchProductArray() {
        for (int gene = 0; gene < GetSize(); gene++) { delete array[gene]; }
    }

    //! return total number of entries (number of genes)
    int GetSize() const { return siteval.GetSize(); }

    //! return total number of branches of the underlying tree
    int GetNbranch() const { return array[0]->GetNbranch(); }

    //! const access to the BranchProduct for the given gene
    const BranchProduct &GetVal(int gene) const {
        // array[gene]->SetMulVal(siteval.GetVal(gene));
        return *array[gene];
    }

    //! non-const access to the BranchProduct for the given gene
    BranchProduct &operator[](int gene) { return *array[gene]; }

    //! global update of the arrays (needs to be called each time the array has
    //! changed)
    void Update() {
        for (int gene = 0; gene < GetSize(); gene++) {
            array[gene]->SetMulVal(siteval.GetVal(gene));
        }
    }

  private:
    const BranchSelector<double> &branchval;
    const Selector<double> &siteval;
    std::vector<BranchProduct *> array;
};


/**
 * \brief The product of two BranchSelector<T>
 *
 * Takes two arguments: a BranchArray l_j and a second BranchArray r_j, and
 * returns, for branch j, the product r_j * l_j.
 *
 * Used in DatedOemgaModel: the branch length (branch j), is
 * equal to the product of Chronogram and BranchRate (time * rate)
 */
template <class T>
class BranchwiseProduct : public SimpleBranchArray<T> {
  public:
    BranchwiseProduct(const BranchSelector<T> &infactor1, const BranchSelector<T> &infactor2)
        : SimpleBranchArray<T>(infactor1.GetTree()), factor1(infactor1), factor2(infactor2) {
        assert(factor1.GetNbranch() == factor2.GetNbranch());
        this->Update();
    }

    //! global update of the branch array
    void Update() {
        for (Tree::BranchIndex branch = 0; branch < this->GetNbranch(); branch++) {
            UpdateBranch(branch);
        }
    }

    //! local update (around a node) of the branch array
    //! Update the branch upstream (parent) and all branches downstream (children)
    void UpdateLocal(Tree::NodeIndex node) {
        // update all branch lengths around this node

        // for all children
        for (auto &child_node : this->GetTree().children(node)) {
            UpdateBranch(this->GetTree().branch_index(child_node));
        }

        if (!this->GetTree().is_root(node)) {
            // for the branch attached to the node
            UpdateBranch(this->GetTree().branch_index(node));
        }
    }

    //! branch update (at a specific branch) of the branch array
    void UpdateBranch(Tree::BranchIndex branch) {
        (*this)[branch] = factor1.GetVal(branch) * factor2.GetVal(branch);
    }

  private:
    const BranchSelector<T> &factor1;
    const BranchSelector<T> &factor2;
};