#pragma once

#include "BranchArray.hpp"
#include "Chronogram.hpp"
#include "NodeArray.hpp"
#include "Random.hpp"
#include "TaxonSet.hpp"
#include "TaxonTraits.hpp"

static int dim_pop_size = 0;
static int dim_mut_rate = 1;
static int dim_gen_time = 2;

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

    //! dimension
    int GetDimensions() const { return dimensions; }

    //! get log prob for a given node
    double GetLogProb(Tree::NodeIndex node) const;

    //! get contrast
    EVector GetContrast(Tree::NodeIndex node) const;

    //! get local log prob for a given node
    double GetLocalLogProb(Tree::NodeIndex node) const;

    //! get global log prob for all nodes
    double GetLogProb() const;

    void ClampLeaves(TaxonTraits const &taxon_traits, TaxonMap const &taxon_map) {
        for (int taxon{0}; taxon < taxon_map.GetNtaxa(); taxon++) {
            Tree::NodeIndex node = taxon_map.TaxonToNode(taxon);
            if (taxon_traits.GenTimePresence() and taxon_traits.GenTimePresence(taxon)) {
                (*this)[node](dim_gen_time) = taxon_traits.GenTime(taxon);
            }
            for (int trait_dim = 0; trait_dim < taxon_traits.GetDim(); trait_dim++) {
                int dim = taxon_traits.TraitDimToMultivariateDim(trait_dim);
                if (taxon_traits.DataPresence(taxon, trait_dim)) {
                    (*this)[node](dim) = taxon_traits.Data(taxon, trait_dim);
                }
            }
        }
        for (Tree::NodeIndex node : GetTree().leaves_root_to_iter()) {
            if (!GetTree().is_leaf(node)) {
                for (int trait_dim = 0; trait_dim < taxon_traits.GetDim(); trait_dim++) {
                    int dim = taxon_traits.TraitDimToMultivariateDim(trait_dim);
                    double sum_children = 0;
                    double sum_branch = 0;
                    for (Tree::NodeIndex child : GetTree().children(node)) {
                        double branch = chronogram.GetVal(tree.branch_index(node));
                        sum_children += (*this)[child](dim) * branch;
                        sum_branch += branch;
                    }
                    (*this)[node](dim) = sum_children / sum_branch;
                }
            }
        }
    }

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

    const Tree &GetTree() const;

    double GetVal(Tree::NodeIndex node) const;
    double &operator[](Tree::NodeIndex node);

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
    explicit BranchProcess(const NodeProcess &innodeprocess, bool arithmetic = false);

    //! global update of the branch array
    void Update();

    //! local update (around a node) of the branch array
    //! Update the branch upstream (parent) and all branches downstream (children)
    void UpdateLocal(Tree::NodeIndex node);

    //! branch update (at a specific branch) of the branch array
    void UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node);

    const NodeProcess &nodeprocess;
    bool arithmetic;
};