#pragma once

#include <algorithm>
#include <cassert>
#include <queue>
#include <vector>
#include "interface.hpp"
#include "nhx-parser.hpp"

// a tree with both a vector of parents and a vector of children
class DoubleVectorTree : public Tree {
    using Tree::NodeIndex;
    std::vector<NodeIndex> parent_;
    std::vector<std::set<NodeIndex>> children_;
    std::vector<std::string> name_;
    NodeIndex root_{0};
    std::vector<NodeIndex> breadth_first_iter;
    std::vector<NodeIndex> breadth_first_inv_iter;

  public:
    explicit DoubleVectorTree(const AnnotatedTree& input_tree) {
        root_ = input_tree.root();
        for (std::size_t i = 0; i < input_tree.nb_nodes(); i++) {
            parent_.push_back(input_tree.parent(i));
            children_.emplace_back(input_tree.children(i).begin(), input_tree.children(i).end());
            name_.push_back(input_tree.tag(i, "name"));
        }
        breadth_first_iter.reserve(nb_nodes());
        breadth_first_inv_iter.resize(nb_nodes());
        std::queue<NodeIndex> bfs_queue;
        bfs_queue.push(root());
        while (!bfs_queue.empty()) {
            NodeIndex node = bfs_queue.front();
            bfs_queue.pop();
            breadth_first_iter.push_back(node);
            for (NodeIndex child : children(node)) { bfs_queue.push(child); }
        }
        assert(breadth_first_iter.size() == nb_nodes());
        std::reverse_copy(
            breadth_first_iter.begin(), breadth_first_iter.end(), breadth_first_inv_iter.begin());
        assert(breadth_first_inv_iter.size() == nb_nodes());
    }

    const std::set<NodeIndex>& children(NodeIndex node) const final { return children_.at(node); }
    NodeIndex parent(NodeIndex node) const final { return parent_.at(node); }
    std::string node_name(NodeIndex node) const final { return name_[node]; }
    NodeIndex root() const final { return root_; }
    std::size_t nb_nodes() const final { return parent_.size(); }
    bool is_root(NodeIndex i) const final { return i == root_; }
    bool is_leaf(NodeIndex i) const final { return children_.at(i).empty(); }
    int nb_branches() const final { return nb_nodes() - 1; }
    BranchIndex branch_index(NodeIndex i) const final { return i - 1; }
    NodeIndex node_index(BranchIndex i) const final { return i + 1; }
    const std::vector<NodeIndex>& root_to_leaves_iter() const final { return breadth_first_iter; }
    const std::vector<NodeIndex>& leaves_root_to_iter() const final { return breadth_first_inv_iter; }
};

std::vector<int> taxa_index_from_parser(TreeParser& parser, const std::vector<std::string>& taxa);
std::unique_ptr<const Tree> make_from_parser(TreeParser& parser);

template <class Element>
std::vector<Element> node_container_from_parser(
    TreeParser& parser, Element (*init)(AnnotatedTree::NodeIndex, const AnnotatedTree&)) {
    using NodeIndex = AnnotatedTree::NodeIndex;
    auto& tree = parser.get_tree();
    std::vector<Element> result;
    for (NodeIndex i = 0; i < NodeIndex(tree.nb_nodes()); i++) { result.push_back(init(i, tree)); }
    return result;
}

template <class Element>
std::vector<Element> branch_container_from_parser(
    TreeParser& parser, Element (*init)(AnnotatedTree::NodeIndex, const AnnotatedTree&)) {
    using NodeIndex = AnnotatedTree::NodeIndex;
    auto& tree = parser.get_tree();
    std::vector<Element> result;
    for (NodeIndex i = 0; i < NodeIndex(tree.nb_nodes()); i++) {
        if (i != tree.root()) {
            result.push_back(init(i, tree));
        } else {
            result.emplace_back();  // default-constructed element for root
        }
    }
    return result;
}

template <class Element>
class TreeElementVector {
    std::vector<Element> v_;
    std::vector<int> index_;

    template <class E>
    friend TreeElementVector<E> taxa_container_from_parser(TreeParser& parser,
        const std::vector<std::string>& taxa,
        E (*init)(AnnotatedTree::NodeIndex, const AnnotatedTree&));

  public:
    using NodeIndex = AnnotatedTree::NodeIndex;

    std::vector<Element>& vector() { return v_; }
    const std::vector<Element>& vector() const { return v_; }
    NodeIndex topology_index(int vector_index) const { return index_.at(vector_index); }
};

template <class Element>
TreeElementVector<Element> taxa_container_from_parser(TreeParser& parser,
    const std::vector<std::string>& taxa,
    Element (*init)(AnnotatedTree::NodeIndex, const AnnotatedTree&)) {
    TreeElementVector<Element> result;
    result.index_ = taxa_index_from_parser(parser, taxa);
    using NodeIndex = AnnotatedTree::NodeIndex;
    auto& tree = parser.get_tree();
    for (size_t i = 0; i < taxa.size(); i++) {
        result.v_.push_back(init(NodeIndex(result.index_.at(i)), tree));
    }
    return result;
}
