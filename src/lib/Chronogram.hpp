#pragma once

#include "BranchArray.hpp"
#include "NodeArray.hpp"

/**
 * \brief A NodeAges
 *
 * Takes one argument: A Tree.
 *
 * Used in DatedNodeMutSelModel: Give an age to each node of the tree, with the constrain that each
 * node is older than its children (may be more than one child) and younger than its parent (at most
 * one parent). The leaves (the youngest) are set to 0, and the root (the oldest) is set to 1.
 */
class NodeAges : public SimpleNodeArray<double> {
  public:
    explicit NodeAges(const Tree& intree, const string& fossils);

    //! find the longest path from this node to the farthest leaf.
    double EccentricityRecursive(Tree::NodeIndex node) const;

    //! check that the constrained are respected
    bool Check() const;

    //! MH move on a node age, respecting the constrains.
    void SlidingMove(Tree::NodeIndex node, double scale);

    bool Unclamped() const { return !clamped; };

  private:
    std::set<Tree::NodeIndex> node_clamped_set{};
    std::unordered_map<Tree::NodeIndex, double> clamped_ages{};
    std::unordered_map<Tree::NodeIndex, double> clamped_lower_bound{};
    std::unordered_map<Tree::NodeIndex, double> clamped_upper_bound{};
    bool clamped{false};
};


/**
 * \brief The Chronogram associated to the underlying NodeAges
 *
 * Takes one argument: a NodeAges.
 *
 * Used in DatedNodeMutSelModel: The time of the branch j is given by the ages of the nodes at the
 * tip of this branch (j).
 */
class Chronogram : public SimpleBranchArray<double> {
  public:
    explicit Chronogram(const NodeAges& innodeages);

    //! global update of the branch array
    void Update();

    //! local update (around a node) of the branch array
    //! Update the branch upstream (parent) and all branches downstream (children)
    void UpdateLocal(Tree::NodeIndex node);

    //! branch update (node and its parent) of the branch array
    void UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node);

  private:
    const NodeAges& nodeages;
};