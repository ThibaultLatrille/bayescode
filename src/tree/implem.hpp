#pragma once

#include <cassert>
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

  public:
    DoubleVectorTree(const AnnotatedTree& input_tree) {
        root_ = input_tree.root();
        for (std::size_t i = 0; i < input_tree.nb_nodes(); i++) {
            parent_.push_back(input_tree.parent(i));
            children_.emplace_back(input_tree.children(i).begin(), input_tree.children(i).end());
            name_.push_back(input_tree.tag(i, "name"));
        }
    }

    const std::set<NodeIndex>& children(NodeIndex node) const final { return children_.at(node); }
    NodeIndex parent(NodeIndex node) const final { return parent_.at(node); }
    std::string node_name(NodeIndex node) const final { return name_[node]; }
    NodeIndex root() const final { return root_; }
    std::size_t nb_nodes() const final { return parent_.size(); }
    bool is_root(NodeIndex i) const final { return i == root_; }
    bool is_leaf(NodeIndex i) const final { return children_.at(i).size() == 0; }
    int nb_branches() const final { return nb_nodes() - 1; }
    BranchIndex branch_index(NodeIndex i) const final { return i - 1; }
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
