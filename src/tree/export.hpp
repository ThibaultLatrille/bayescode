#pragma once

#include <unordered_map>
#include <vector>
#include "interface.hpp"

// a tree with both a vector of parents and a vector of children
class ExportTree {
    using TagName = std::string;
    using TagValue = std::string;
    using Node = std::unordered_map<TagName, TagValue>;

    std::vector<Node> nodes_;
    Tree const &tree;

  public:
    explicit ExportTree(Tree const &intree) : tree{intree} {
        nodes_.resize(tree.nb_nodes());
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
            set_tag(node, "name", tree.node_name(node));
        }
    }

    Tree const &GetTree() { return tree; };

    void set_tag(Tree::NodeIndex node, TagName const &tag, TagValue const &value) {
        nodes_.at(node)[tag] = value;
    }

    TagValue tag(Tree::NodeIndex node, TagName const &tag) const {
        if (nodes_.at(node).count(tag) != 0) {
            return nodes_.at(node).at(tag);
        } else {
            return "";
        }
    }

    std::string recursive_string(Tree::NodeIndex node) const {
        std::string newick;

        if (not tree.children(node).empty()) {
            // It's an internal node
            newick += "(";
            for (auto const child : tree.children(node)) {
                newick += recursive_string(child) + ",";
            }
            newick.pop_back();
            newick += ")";
        }
        auto node_annotation = nodes_.at(node);
        if (node_annotation.count("name") != 0) {
            newick += nodes_.at(node).at("name");
            node_annotation.erase("name");
        }
        if (node_annotation.count("length") != 0) {
            newick += ":" + nodes_.at(node).at("length");
            node_annotation.erase("length");
        }
        if (not node_annotation.empty()) {
            newick += "[&&NHX";
            for (auto &it : node_annotation) { newick += ":" + it.first + "=" + it.second; }
            newick += "]";
        }
        return newick;
    }

    std::string as_string() const { return recursive_string(tree.root()) + "; "; }
};
