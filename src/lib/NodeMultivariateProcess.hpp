#pragma once

#include <cassert>
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
        const Chronogram &inchrono, const EMatrix &inprecision_matrix, const EVector &inroot_mean)
        : SimpleNodeArray<EVector>(inchrono.GetTree()),
          chronogram(inchrono),
          dimensions(inroot_mean.size()),
          precision_matrix(inprecision_matrix),
          root_mean(inroot_mean) {
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
            (*this)[node] = EVector::Zero(dimensions);
        }
        Sample();
    }

    //! sample all entries from prior
    void Sample() {
        (*this)[GetTree().root()].Zero(dimensions);
        for (int dim = 0; dim < dimensions; dim++) {
            (*this)[GetTree().root()](dim) = root_mean(dim) + GetSigma(dim) * Random::sNormal();
        };
        SampleRecursive(GetTree().root());
    }

    //! sample recursively a node from prior
    void SampleRecursive(Tree::NodeIndex node) {
        for (auto const &child : GetTree().children(node)) {
            for (int dim = 0; dim < dimensions; dim++) {
                (*this)[child](dim) =
                    this->GetVal(node)(dim) + GetSigmaScaled(child, dim) * Random::sNormal();
            };
            SampleRecursive(child);
        }
    }

    //! variance of the pro recursively a node from prior
    double GetSigma(int dimension) const {
        return sqrt(1.0 / precision_matrix(dimension, dimension));
    };

    double GetSigmaScaled(Tree::NodeIndex node, int dimension) const {
        return GetSigma(dimension) *
               sqrt(chronogram.GetVal(chronogram.GetTree().branch_index(node)));
    };

    //! dimension
    int GetDimensions() const { return dimensions; };

    //! get log prob for a given node
    double GetLogProb(Tree::NodeIndex node) const {
        return Random::logNormalDensity(GetContrast(node), precision_matrix);
    }

    //! get contrast
    EVector GetContrast(Tree::NodeIndex node) const {
        if (GetTree().is_root(node)) {
            return this->GetVal(node) - root_mean;
        } else {
            return (this->GetVal(node) - this->GetVal(GetTree().parent(node))) /
                   sqrt(chronogram.GetVal(chronogram.GetTree().branch_index(node)));
        }
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
    const EVector &root_mean;
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