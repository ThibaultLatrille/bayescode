#pragma once

#include <set>
#include <string>

class Tree {
  public:
    using NodeIndex = int;  // guaranteed to be [0... nb_nodes() - 1] U {no_parent}
    const NodeIndex no_parent{-1};
    using BranchIndex = int;  // guaranteed to be [0... nb_nodes() - 2]

    // topology-related methods
    virtual const std::set<NodeIndex>& children(NodeIndex) const = 0;
    virtual NodeIndex parent(NodeIndex) const = 0;  // parent(root()) returns no_parent
    virtual std::string node_name(NodeIndex) const = 0;
    virtual NodeIndex root() const = 0;
    virtual std::size_t nb_nodes() const = 0;
    virtual bool is_leaf(NodeIndex) const = 0;  // (for convenience)
    virtual bool is_root(NodeIndex) const = 0;  // (for convenience)
    // virtual const std::set<NodeIndex>& neighbors(NodeIndex) const = 0;  // for unrooted searches,
    // non-negligible cost
    virtual int nb_branches() const = 0;
    virtual BranchIndex branch_index(NodeIndex) const = 0;
    virtual NodeIndex node_index(BranchIndex) const = 0;

    virtual ~Tree() = default;
};
