#pragma once

#include <cassert>
#include "BranchArray.hpp"
#include "NodeArray.hpp"

/**
 * \brief A NodeAges
 *
 * Takes one argument: A Tree.
 *
 * Used in DatedOmegaModel: Give an age to each node of the tree, with the constrain that each node
 * is older than its children (may be more than one child) and younger than its parent (at most one
 * parent). The leaves (the youngest) are set to 0, and the root (the oldest) is set to 1.
 */
class NodeAges : public SimpleNodeArray<double> {
  public:
    explicit NodeAges(const Tree& intree) : SimpleNodeArray<double>(intree) { UniformSample(); }

    //! sample all entries uniformly by finding the longest path to the leaves.
    void UniformSample() {
        double eccentricity = EccentricityRecursive(GetTree().root());

        for (Tree::NodeIndex node{0}; node < GetNnode(); node++) {
            // renormalize ages so that root age is 1 by definition
            (*this)[node] /= eccentricity;
        }
        Check();
    }

    //! find the longest path from this node to the farthest leaf.
    double EccentricityRecursive(Tree::NodeIndex node) {
        if (GetTree().is_leaf(node)) {
            // set leaf nodes at age 0
            (*this)[node] = 0.0;
        } else {
            // proceed from leaves to root using recursive algorithm
            double max_eccent = 0.0;
            for (auto const& child : GetTree().children(node)) {
                double eccent = EccentricityRecursive(child) + 1;
                if (eccent > max_eccent) { max_eccent = eccent; }
            }

            (*this)[node] = max_eccent;
        }
        return GetVal(node);
    };

    //! check that the constrained are respected
    void Check() const {
        for (Tree::NodeIndex node{0}; node < GetNnode(); node++) {
            if (GetTree().is_root(node)) {
                assert(this->GetVal(node) == 1.0);
            } else if (GetTree().is_leaf(node)) {
                assert(this->GetVal(node) == 0.0);
            } else {
                assert(this->GetVal(GetTree().parent(node)) >= this->GetVal(node));
            }
        }
    }

    //! MH move on a node age, respecting the constrains.
    void SlidingMove(Tree::NodeIndex node, double scale) {
        assert(!GetTree().is_root(node));
        assert(!GetTree().is_leaf(node));

        double lower_bound = 0.0;
        for (auto const& child : GetTree().children(node)) {
            if (this->GetVal(child) > lower_bound) { lower_bound = this->GetVal(child); }
        }
        double upper_bound = this->GetVal(GetTree().parent(node));

        assert(upper_bound >= lower_bound);

        double x = this->GetVal(node) + scale * (upper_bound - lower_bound);

        while ((x < lower_bound) || (x > upper_bound)) {
            if (x < lower_bound) { x = 2 * lower_bound - x; }
            if (x > upper_bound) { x = 2 * upper_bound - x; }
        }

        (*this)[node] = x;
    }
};

/**
 * \brief The Chronogram associated to the underlying NodeAges
 *
 * Takes one argument: a NodeAges.
 *
 * Used in DatedOemgaModel: The time of the branch j is given by the ages of the nodes at the tip
 * of this branch (j).
 */
class Chronogram : public SimpleBranchArray<double> {
  public:
    explicit Chronogram(const NodeAges& innodeages)
        : SimpleBranchArray<double>(innodeages.GetTree()), nodeages(innodeages) {
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
        for (auto const& child : this->GetTree().children(node)) { UpdateBranch(node, child); }

        // for the branch attached to the node
        if (!this->GetTree().is_root(node)) { UpdateBranch(this->GetTree().parent(node), node); }
    }

    //! branch update (node and its parent) of the branch array
    void UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node) {
        (*this)[this->GetTree().branch_index(node)] =
            nodeages.GetVal(parent) - nodeages.GetVal(node);
    }

  private:
    const NodeAges& nodeages;
};