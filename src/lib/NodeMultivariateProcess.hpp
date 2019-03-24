#pragma once

#include <cassert>
#include <cmath>
#include "BranchArray.hpp"
#include "Chronogram.hpp"
#include "NodeArray.hpp"
#include "Random.hpp"

/**
 * \brief A brownian univariate NodeProcess
 *
 * Takes three arguments: A Chronogram (time associated to branches), a sigma (variance of the
 * process) and the expectation of the process at the root.
 *
 * Used in DatedOmegaModel: each node is dependent on its parent.
 * The value at the node is normally distributed of mean equal to the parent node, and the variance
 * is given by the time of the branch (Chronogram) multiplied by sigma^2.
 */
class NodeMultivariateProcess : public SimpleNodeArray<EVector> {
  public:
    NodeMultivariateProcess(
        const Chronogram &inchrono, const EMatrix &inprecision_matrix, int indimensions)
        : SimpleNodeArray<EVector>(inchrono.GetTree()),
          chronogram(inchrono),
          dimensions(indimensions),
          precision_matrix(inprecision_matrix) {
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
            (*this)[node] = EVector::Zero(dimensions);
        }
    }

    //! variance of the pro recursively a node from prior
    double GetSigma(int dimension) const {
        return sqrt(1.0 / precision_matrix(dimension, dimension));
    };

    //! dimension
    int GetDimensions() const { return dimensions; };

    //! get log prob for a given node
    double GetLogProb(Tree::NodeIndex node) const {
        if (GetTree().is_root(node)) {
            return 0.0;
        } else {
            return Random::logNormalDensity(GetContrast(node), precision_matrix);
        }
    }

    //! get contrast
    EVector GetContrast(Tree::NodeIndex node) const {
        assert(!GetTree().is_root(node));
        return (this->GetVal(node) - this->GetVal(GetTree().parent(node))) /
               sqrt(chronogram.GetVal(GetTree().branch_index(node)));
    }

    //! get local log prob for a given node
    double GetLocalLogProb(Tree::NodeIndex node) const {
        double tot = GetLogProb(node);
        for (auto const &child : GetTree().children(node)) { tot += GetLogProb(child); }
        return tot;
    }

    //! get global log prob for all nodes
    double GetLogProb() const {
        double tot = 0;
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
            if (!GetTree().is_root(node)) { tot += GetLogProb(node); }
        }
        return tot;
    }

  protected:
    const Chronogram &chronogram;
    int dimensions;
    const EMatrix &precision_matrix;
};

class NodeProcess {
  public:
    NodeProcess(NodeMultivariateProcess &innode_multivariate, int indimension)
        : dimension(indimension), node_multivariate(innode_multivariate) {}

    //! variance of the pro recursively a node from prior
    double GetSigma() const { return node_multivariate.GetSigma(dimension); };

    const Tree &GetTree() const { return node_multivariate.GetTree(); }

    double GetVal(Tree::NodeIndex node) const { return node_multivariate.GetVal(node)(dimension); }
    double &operator[](Tree::NodeIndex node) { return node_multivariate[node](dimension); }

    double GetExpVal(Tree::NodeIndex node) const { return exp(GetVal(node)); }

    void SlidingMove(Tree::NodeIndex node, double m) { node_multivariate[node](dimension) += m; }

    void SlidingMove(double m) {
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
            SlidingMove(node, m);
        }
    }

  protected:
    int dimension;
    NodeMultivariateProcess &node_multivariate;
};

/**
 * \brief The branch process associated to the underlying NodeProcess
 *
 * Takes one argument: a NodeProcess.
 *
 * Used in DatedOemgaModel: The value of the process at branch j is the exponential of the half-sum
 * of the NodeProcess, for the nodes at the tip of branch j.
 */
class BranchProcess : public SimpleBranchArray<double> {
  public:
    //! Constructor (with only the tree given as argument)
    explicit BranchProcess(const NodeProcess &innodeprocess)
        : SimpleBranchArray<double>(innodeprocess.GetTree()), nodeprocess{innodeprocess} {
        Update();
    }

    //! global update of the branch array
    void Update() {
        for (Tree::NodeIndex node{0}; node < Tree::NodeIndex(this->GetTree().nb_nodes()); node++) {
            if (!this->GetTree().is_root(node)) {
                this->UpdateBranch(this->GetTree().parent(node), node);
            }
        }
    }

    //! local update (around a node) of the branch array
    //! Update the branch upstream (parent) and all branches downstream (children)
    void UpdateLocal(Tree::NodeIndex node) {
        // update all branch lengths around this node

        // for all children
        for (auto const &child : this->GetTree().children(node)) { UpdateBranch(node, child); }

        // for the branch attached to the node
        if (!this->GetTree().is_root(node)) { UpdateBranch(this->GetTree().parent(node), node); }
    }

    //! branch update (at a specific branch) of the branch array
    void UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node) {
        (*this)[this->GetTree().branch_index(node)] =
            (nodeprocess.GetExpVal(parent) + nodeprocess.GetExpVal(node)) / 2;
        // For geodesic, use :
        // = (nodeprocess.GetExpVal(parent) - nodeprocess.GetExpVal(node)) /
        // (nodeprocess->GetVal(parent) - nodeprocess->GetVal(node));
    }

    const NodeProcess &nodeprocess;
};

class BranchWiseMultivariateProcess : public SimpleBranchArray<EVector> {
  public:
    BranchWiseMultivariateProcess(
        const Chronogram &inchrono, const EMatrix &inprecision_matrix, int indimensions)
        : SimpleBranchArray<EVector>(inchrono.GetTree()),
          chronogram(inchrono),
          dimensions(indimensions),
          precision_matrix(inprecision_matrix) {
        root_process = EVector::Zero(dimensions);
        for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
            (*this)[branch] = EVector::Zero(dimensions);
        }
    }

    //! dimension
    int GetDimensions() const { return dimensions; };

    //! variance of the pro recursively a node from prior
    double GetSigma(int dimension) const {
        return sqrt(1.0 / precision_matrix(dimension, dimension));
    };

    //! get log prob for a given node
    double GetLogProb(Tree::BranchIndex branch) const {
        return Random::logNormalDensity(GetContrast(branch), precision_matrix);
    }

    //! get contrast
    EVector GetContrast(Tree::BranchIndex branch) const {
        double distance = chronogram.GetVal(branch) / 2;
        Tree::NodeIndex node = GetTree().node_index(branch);
        assert(!GetTree().is_root(node));
        Tree::NodeIndex parent_node = GetTree().parent(node);
        if (GetTree().is_root(parent_node)) {
            return (root_process - this->GetVal(branch)) / sqrt(distance);
        } else {
            Tree::BranchIndex parent_branch = GetTree().branch_index(parent_node);
            distance += chronogram.GetVal(parent_branch) / 2;
            return (this->GetVal(parent_branch) - this->GetVal(branch)) / sqrt(distance);
        }
    }

    //! get local log prob for a given node
    double GetLocalLogProb(Tree::NodeIndex node) const {
        double tot = GetLogProb(GetTree().branch_index(node));
        for (Tree::NodeIndex const &child : GetTree().children(node)) {
            tot += GetLogProb(GetTree().branch_index(child));
        }
        return tot;
    }

    //! get global log prob for all branches
    double GetLogProb() const {
        double tot = 0;
        for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
            tot += GetLogProb(branch);
        }
        return tot;
    }

    EVector root_process;
    const Chronogram &chronogram;
    int dimensions;
    const EMatrix &precision_matrix;
};

class LeafMultivariateProcess : public SimpleNodeArray<EVector> {
  public:
    LeafMultivariateProcess(
        const BranchWiseMultivariateProcess &inbranchwiseprocess, TaxonSet const &taxon)
        : SimpleNodeArray<EVector>(inbranchwiseprocess.GetTree()),
          branchwiseprocess(inbranchwiseprocess) {
        index_table = taxon.get_index_table(&branchwiseprocess.GetTree());
        for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
            (*this)[branch] = EVector::Zero(branchwiseprocess.GetDimensions());
        }
    }

    //! get log prob for a given taxon
    double GetTaxonLogProb(int taxon) const { return GetLogProb(index_table[taxon]); }

    //! get log prob for a given node
    double GetLogProb(Tree::NodeIndex node) const {
        assert(GetTree().is_leaf(node));
        return Random::logNormalDensity(GetContrast(node), branchwiseprocess.precision_matrix);
    }

    //! get contrast
    EVector GetContrast(Tree::NodeIndex node) const {
        Tree::BranchIndex parent_branch = GetTree().branch_index(node);
        double distance = branchwiseprocess.chronogram.GetVal(parent_branch) / 2;
        return (branchwiseprocess.GetVal(parent_branch) - this->GetVal(node)) / sqrt(distance);
    }

    double GetTheta(Tree::NodeIndex node) const {
        assert(GetTree().is_leaf(node));
        return exp(this->GetVal(node).sum());
    }

    //! get global log prob for all branches
    double GetLogProb() const {
        double tot = 0;
        for (Tree::NodeIndex node = 0; node < static_cast<Tree::NodeIndex>(GetTree().nb_nodes());
             node++) {
            if (tree.is_leaf(node)) { tot += GetLogProb(node); }
        }
        return tot;
    }

  protected:
    std::vector<int> index_table;
    const BranchWiseMultivariateProcess &branchwiseprocess;
};

class BranchWiseProcess : public SimpleBranchArray<double> {
  public:
    BranchWiseProcess(BranchWiseMultivariateProcess &inbranchwise_multivariate, int indimension)
        : SimpleBranchArray<double>(inbranchwise_multivariate.GetTree()),
          dimension(indimension),
          branchwise_multivariate(inbranchwise_multivariate) {}

    //! variance of the pro recursively a node from prior
    double GetSigma() const { return branchwise_multivariate.GetSigma(dimension); };

    void SlidingMove(Tree::BranchIndex branch, double m) {
        branchwise_multivariate[branch](dimension) += m;
        this->array[branch] = exp(branchwise_multivariate[branch](dimension));
    }

    double GetRootVal() { return exp(branchwise_multivariate.root_process(dimension)); }

    void SlidingRootMove(double m) { branchwise_multivariate.root_process(dimension) += m; }

    void Update() {
        std::vector<double> array(GetTree().nb_branches(), 0.0);
        for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
            this->array[branch] = exp(branchwise_multivariate[branch](dimension));
        }
    }

  protected:
    int dimension;
    BranchWiseMultivariateProcess &branchwise_multivariate;
};
