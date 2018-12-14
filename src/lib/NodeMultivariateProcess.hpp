#pragma once

#include <cassert>
#include "BranchArray.hpp"
#include "Chronogram.hpp"
#include "NodeArray.hpp"
#include "Random.hpp"

class CovMatrix : public EMatrix {
  public:
    explicit CovMatrix(int indimensions)
        : EMatrix(indimensions, indimensions), dimensions(indimensions) {
        for (int row = 1; row < dimensions; row++) {
            for (int col = 1; col < dimensions; col++) {
                (*this)(row, col) = 0.0;
            };
            (*this)(row, row) = 1.0;
        };
    }

    int GetDimensions() const {
        return dimensions;
    }

    void SlidingMove(int i, int j, double v) {
        if (i == j) {
            (*this)(i, j) += v;
            // easy to assert that sigma
        }
    };

  protected:
    int dimensions;
};


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
    NodeMultivariateProcess(const Chronogram &inchrono, const CovMatrix &incov_matrix,
        const EVector &inroot_mean)
        : SimpleNodeArray<EVector>(inchrono.GetTree()),
          chronogram(inchrono),
          dimensions(incov_matrix.GetDimensions()),
          cov_matrix(incov_matrix),
          root_mean(inroot_mean) {
        assert(root_mean.size() == dimensions);
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
            (*this)[node] = EVector::Zero(dimensions);
        }
        Sample();
    }

    //! sample all entries from prior
    void Sample() {
        (*this)[GetTree().root()].Zero(dimensions);
        for (int dim = 1; dim < dimensions; dim++) {
            double test = root_mean(dim);
            double cov =  cov_matrix(dim, dim);
            (*this)[GetTree().root()](dim) =
                test * exp(cov * Random::sNormal());
        };
        SampleRecursive(GetTree().root());
    }

    //! sample recursively a node from prior
    void SampleRecursive(Tree::NodeIndex node) {
        for (auto const &child : GetTree().children(node)) {
            for (int dim = 0; dim < dimensions; dim++) {
                (*this)[child](dim) =
                    this->GetVal(node)(dim) * exp(cov_matrix(dim, dim) * Random::sNormal());
            };
            SampleRecursive(child);
        }
    }

    //! variance of the pro recursively a node from prior
    double GetSigma(int dimension) const { return sqrt(cov_matrix(dimension, dimension)); };

    //! dimension
    int GetDimensions() const { return dimensions; };

    //! get log prob for a given node
    double GetLogProb(Tree::NodeIndex node) const {
        return Random::logNormalDensity(GetContrast(node), cov_matrix);
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
            if (!GetTree().is_root(node)){
                tot += GetLogProb(node);
            }
        }
        return tot;
    }

  protected:
    const Chronogram &chronogram;
    int dimensions;
    const CovMatrix &cov_matrix;
    const EVector &root_mean;
};

class NodeProcess {
  public:
    NodeProcess(NodeMultivariateProcess &innode_multivariate, int indimension)
        : dimension(indimension), node_multivariate(innode_multivariate)  {}

    //! variance of the pro recursively a node from prior
    double GetSigma() const { return node_multivariate.GetSigma(dimension); };

    const Tree &GetTree() const {
        return node_multivariate.GetTree();
    }

    double GetVal(Tree::NodeIndex node) const {
        return node_multivariate.GetVal(node)(dimension);
    }

    void SlidingMove(Tree::NodeIndex node, double m) {
        node_multivariate[node](dimension) += m;
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
            (exp(nodeprocess.GetVal(parent)) + exp(nodeprocess.GetVal(node))) / 2;
        // For geodesic, use :
        // = (exp(nodeprocess->GetVal(parent)) - exp(nodeprocess->GetVal(node))) /
        // (nodeprocess->GetVal(parent) - nodeprocess->GetVal(node));
    }

    //! get sum over the branch array
    double GetSum() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) { m1 += GetVal(i); }
        return m1;
    }

    //! get mean over the branch array
    double GetMean() const { return GetSum() / GetNbranch(); }

    //! get variance over the branch array
    double GetVar() const {
        double m1 = GetMean();
        double m2 = 0;
        for (int i = 0; i < GetNbranch(); i++) { m2 += GetVal(i) * GetVal(i); }
        m2 /= GetNbranch();
        m2 -= m1 * m1;
        return m2;
    }
    const NodeProcess &nodeprocess;
};