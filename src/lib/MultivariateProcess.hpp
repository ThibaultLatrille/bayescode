#pragma once

#include "BranchArray.hpp"
#include "Chronogram.hpp"
#include "NodeArray.hpp"
#include "Random.hpp"
#include "TaxonSet.hpp"

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
        const Chronogram &inchrono, const EMatrix &inprecision_matrix, int indimensions);

    //! variance of the pro recursively a node from prior
    double GetSigma(int dimension) const;
    ;

    //! dimension
    int GetDimensions() const { return dimensions; };

    //! get log prob for a given node
    double GetLogProb(Tree::NodeIndex node) const;

    //! get contrast
    EVector GetContrast(Tree::NodeIndex node) const;

    //! get local log prob for a given node
    double GetLocalLogProb(Tree::NodeIndex node) const;

    //! get global log prob for all nodes
    double GetLogProb() const;

  protected:
    const Chronogram &chronogram;
    int dimensions;
    const EMatrix &precision_matrix;
};

class NodeProcess {
  public:
    NodeProcess(NodeMultivariateProcess &innode_multivariate, int indimension);

    //! variance of the pro recursively a node from prior
    double GetSigma() const;
    ;

    const Tree &GetTree() const;

    double GetVal(Tree::NodeIndex node) const;
    double &operator[](Tree::NodeIndex node);

    double GetExpVal(Tree::NodeIndex node) const;

    void SlidingMove(Tree::NodeIndex node, double m);

    void SlidingMove(double m);

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
    explicit BranchProcess(const NodeProcess &innodeprocess);

    //! global update of the branch array
    void Update();

    //! local update (around a node) of the branch array
    //! Update the branch upstream (parent) and all branches downstream (children)
    void UpdateLocal(Tree::NodeIndex node);

    //! branch update (at a specific branch) of the branch array
    void UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node);

    const NodeProcess &nodeprocess;
};

class BranchWiseMultivariateProcess : public SimpleBranchArray<EVector> {
  public:
    BranchWiseMultivariateProcess(
        const Chronogram &inchrono, const EMatrix &inprecision_matrix, int indimensions);

    //! dimension
    int GetDimensions() const { return dimensions; };

    double GetBranchLogProb(Tree::BranchIndex branch) const;

    //! get log prob for a given node
    double GetLocalBranchLogProb(Tree::BranchIndex branch) const;

    //! get contrast
    EVector GetContrast(Tree::BranchIndex branch) const;

    //! get local log prob for a given node
    double GetLocalNodeLogProb(Tree::NodeIndex node) const;

    //! get global log prob for all branches
    double GetLogProb() const;

    EVector root_process;
    const Chronogram &chronogram;
    int dimensions;
    const EMatrix &precision_matrix;
};

class LeafMultivariateProcess : public SimpleNodeArray<EVector> {
  public:
    LeafMultivariateProcess(
        const BranchWiseMultivariateProcess &inbranchwiseprocess, TaxonSet const &taxon);

    //! get log prob for a given taxon
    double GetTaxonLogProb(int taxon) const;

    //! get log prob for a given node
    double GetLogProb(Tree::NodeIndex node) const;

    //! get contrast
    EVector GetContrast(Tree::NodeIndex node) const;

    double GetTheta(Tree::NodeIndex node) const;

    //! get global log prob for all branches
    double GetLogProb() const;

  protected:
    std::vector<int> reverse_index_table;
    const BranchWiseMultivariateProcess &branchwiseprocess;
};

class BranchWiseProcess : public SimpleBranchArray<double> {
  public:
    BranchWiseProcess(BranchWiseMultivariateProcess &inbranchwise_multivariate, int indimension);

    void SlidingMove(Tree::BranchIndex branch, double m);

    double GetRootVal() const { return root_value; }

    void SlidingRootMove(double m);

    void Update();

  protected:
    double root_value;
    int dimension;
    BranchWiseMultivariateProcess &branchwise_multivariate;
};
