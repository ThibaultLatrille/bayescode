/*Copyright or Copr. Centre National de la Recherche Scientifique (CNRS) (2018)
Contributors:
- Vincent Lanore <vincent.lanore@gmail.com>

This software is a computer program whose purpose is to provide a header-only standalone parser for
NHX (New Hampshire Extended) phylogenetic trees.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/

#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

/*
====================================================================================================
  ~*~ Pure interfaces ~*~
==================================================================================================*/
class AnnotatedTree {
  public:
    using NodeIndex = int;
    using TagName = std::string;
    using TagValue = std::string;
    using ChildrenList = std::vector<NodeIndex>;

    virtual const ChildrenList& children(NodeIndex) const = 0;
    virtual NodeIndex parent(NodeIndex) const = 0;
    virtual NodeIndex root() const = 0;
    virtual std::size_t nb_nodes() const = 0;
    virtual TagValue tag(NodeIndex, TagName) const = 0;

    virtual std::string as_string() const = 0;
    virtual std::vector<std::string> descendant_leaves(NodeIndex) const = 0;
    virtual bool operator==(const AnnotatedTree& other) const = 0;

    virtual ~AnnotatedTree() = default;
};

class TreeParser {
  public:
    virtual const AnnotatedTree& get_tree() const = 0;
    virtual ~TreeParser() = default;
};

/*
====================================================================================================
  ~*~ Implementations ~*~
==================================================================================================*/
class DoubleListAnnotatedTree : public AnnotatedTree {
  public:
    using Node = std::unordered_map<std::string, std::string>;

    std::vector<Node> nodes_;

    // invariant: same length as nodes
    // invariant: only one node has parent -1
    // element i is the index of parent of i in nodes (-1 for root)
    std::vector<int> parent_;

    // invariant: consistent with parent
    // element i is a vector of indices corresponding to the children of node i
    std::vector<std::vector<int>> children_;

    // invariant: node with index root is only node with parent -1
    NodeIndex root_;

  public:
    const ChildrenList& children(NodeIndex node) const final { return children_.at(node); }

    NodeIndex parent(NodeIndex node) const final { return parent_.at(node); }

    NodeIndex root() const final { return root_; }

    std::size_t nb_nodes() const final { return nodes_.size(); }

    TagValue tag(NodeIndex node, TagName tag) const final {
        if (nodes_.at(node).count(tag) != 0) {
            return nodes_.at(node).at(tag);
        } else {
            return "";
        }
    }

    std::string recursive_string(NodeIndex node) const {
        std::string newick;

        if (not children(node).empty()) {
            // It's an internal node
            newick += "(";
            for (auto const child : children(node)) { newick += recursive_string(child) + ","; };
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
            for (auto& it : node_annotation) { newick += ":" + it.first + "=" + it.second; }
            newick += "]";
        }
        return newick;
    }

    std::string as_string() const final { return recursive_string(root()) + "; "; }

    std::vector<std::string> descendant_leaves(NodeIndex node) const final {
        std::vector<std::string> leaves(0);
        if (children(node).empty()) {
            leaves.emplace_back(tag(node, "name"));
        } else {
            for (auto const& child : children(node)) {
                auto v2 = descendant_leaves(child);
                if (not v2.empty()) { leaves.insert(leaves.end(), v2.begin(), v2.end()); }
            }
        }
        return leaves;
    }

    int recursive_diff(NodeIndex node, const AnnotatedTree& other, NodeIndex other_node) const {
        int diff = 0;

        if (children(node).size() != other.children(other_node).size()) {
            diff++;
        } else {
            for (auto const& child : children(node)) {
                auto child_leaves = descendant_leaves(child);
                std::sort(child_leaves.begin(), child_leaves.end());
                bool no_counterpart = true;
                for (auto const& other_child : other.children(other_node)) {
                    auto other_leaves = other.descendant_leaves(other_child);
                    std::sort(other_leaves.begin(), other_leaves.end());

                    if (tag(child, "name") == other.tag(other_child, "name") and
                        tag(child, "length") == other.tag(other_child, "length") and
                        child_leaves == other_leaves) {
                        no_counterpart = false;
                        diff += recursive_diff(child, other, other_child);
                    }
                }
                if (no_counterpart) { diff++; }
            }
        }

        for (auto& it : nodes_.at(node)) {
            if (it.second != other.tag(other_node, it.first)) { diff++; };
        }

        return diff;
    }

    bool operator==(const AnnotatedTree& other) const final {
        if (nb_nodes() != other.nb_nodes()) { return false; };
        return recursive_diff(root(), other, other.root()) == 0;
    }
};

/*================================================================================================*/
struct NHXParserException : public std::runtime_error {
    NHXParserException(std::string s = "") : std::runtime_error(s) {}
};

/*================================================================================================*/
class NHXParser : public TreeParser {
    // list of tokens for lexer
    enum TokenType {
        OpenParenthesis,
        CloseParenthesis,
        Colon,
        Semicolon,
        Comma,
        Equal,
        NHXOpen,
        CommentOpen,
        BracketClose,
        Identifier,
        Invalid
    };
    std::map<TokenType, std::string> token_names{{OpenParenthesis, "OpenParenthesis"},
        {CloseParenthesis, "CloseParenthesis"}, {Colon, "Colon"}, {Semicolon, "Semicolon"},
        {Comma, "Comma"}, {Equal, "Equal"}, {NHXOpen, "NHXOpen"}, {CommentOpen, "CommentOpen"},
        {BracketClose, "BracketClose"}, {Identifier, "Identifier"}, {Invalid, "Invalid"}};
    std::map<TokenType, std::regex> token_regexes{{OpenParenthesis, std::regex("\\(")},
        {CloseParenthesis, std::regex("\\)")}, {Colon, std::regex(":")},
        {Semicolon, std::regex(";")}, {Comma, std::regex(",")}, {Equal, std::regex("=")},
        {NHXOpen, std::regex("\\[&&NHX:")}, {CommentOpen, std::regex("\\[")},
        {BracketClose, std::regex("\\]")}, {Identifier, std::regex("[a-zA-Z0-9._-]+")}};
    using Token = std::pair<TokenType, std::string>;  // first: index of token, second: token value

    // input/output
    DoubleListAnnotatedTree tree;

    // state during parsing
    using scit = std::string::const_iterator;
    scit it;
    Token next_token{Invalid, ""};
    std::string input{""};
    int next_node{0};

    [[noreturn]] void error(std::string s) {
        std::stringstream ss;
        ss << s;
        ss << "Error at position " << std::distance(scit(input.begin()), it) << ":\n";
        bool at_begining = it - 15 <= input.begin();
        bool at_end = it + 15 >= input.end();
        ss << "\t" << (at_begining ? "" : "...")
           << std::string(at_begining ? input.begin() : it - 15, at_end ? input.end() : it + 15)
           << (at_end ? "" : "...") << "\n";
        ss << "\t" << (at_begining ? "" : "   ")
           << std::string(at_begining ? std::distance(scit(input.begin()), it) : 15, ' ') << "^\n";
        throw NHXParserException(ss.str());
    }

    std::string expect(TokenType type) {
        find_token();
        if (next_token.first != type) {
            error("Error: expected token " + token_names.at(type) + " but got token " +
                  token_names.at(next_token.first) + "(" + next_token.second + ") instead.\n");
        } else {
            return next_token.second;
        }
    }

    // lexer
    void find_token() {
        while (std::isspace(*it)) { it++; }
        int token_number{0};
        for (auto token_regex : token_regexes) {
            std::smatch m;
            auto r = std::regex_search(it, scit(input.end()), m, token_regex.second);
            if (r and m.prefix() == "") {
                if (token_regex.first == CommentOpen) {  // support of comments
                    std::string comment_close{"]"};
                    it = std::search(
                             it, scit(input.end()), comment_close.begin(), comment_close.end()) +
                         1;
                    find_token();
                    return;
                }
                next_token = Token{token_regex.first, m[0]};
                it += std::string(m[0]).size();
                // std::cout << "found token " << m[0] << std::endl;
                return;
            }
            token_number++;
        }
        next_token = Token{Invalid, it == scit(input.end())
                                        ? "end of input"
                                        : ("token starting with " + std::string(it, it + 1))};
    }

    // parser
    void node_nothing(int number, int parent) {
        tree.nodes_.emplace_back();
        tree.parent_.push_back(parent);
        tree.children_.emplace_back();
        if (parent != -1) { tree.children_.at(parent).push_back(number); }

        find_token();
        switch (next_token.first) {
            case Identifier:
                tree.nodes_[number]["name"] = next_token.second;
                node_name(number, parent);
                break;
            case Colon: node_length(number, parent); break;
            case NHXOpen: data(number, parent); break;
            case OpenParenthesis:
                next_node++;
                node_nothing(next_node, number);
                break;
            default: node_end(parent);
        }
    }

    void node_name(int number, int parent) {
        find_token();
        switch (next_token.first) {
            case Colon: node_length(number, parent); break;
            case NHXOpen: data(number, parent); break;
            case Identifier:
                tree.nodes_[number]["name"] = next_token.second;
                node_name(number, parent);
                break;
            default: node_end(parent);
        }
    }

    void node_length(int number, int parent) {
        tree.nodes_[number]["length"] = expect(Identifier);

        find_token();
        switch (next_token.first) {
            case NHXOpen: data(number, parent); break;
            default: node_end(parent);
        }
    }

    void node_end(int parent) {
        switch (next_token.first) {
            case Comma:
                next_node++;
                node_nothing(next_node, parent);
                break;
            case CloseParenthesis: {
                if (parent != -1) { node_name(parent, tree.parent_.at(parent)); }
                break;
            }
            case Semicolon: break;
            default: error("Error: unexpected " + next_token.second + '\n');
        }
    }

    void data(int number, int parent) {
        find_token();
        if (next_token.first == BracketClose) {
            find_token();
            node_end(parent);
        } else if (next_token.first == Identifier) {
            std::string tag = next_token.second;
            expect(Equal);
            tree.nodes_[number][tag] = expect(Identifier);
            data(number, parent);
        } else if (next_token.first == Colon) {
            data(number, parent);
        } else {
            error("Error: improperly formatted contents in NHX data. Found unexpected " +
                  next_token.second + '\n');
        }
    }

  public:
    NHXParser(std::istream& is) {
        tree = DoubleListAnnotatedTree();
        tree.root_ = 0;
        input = std::string(std::istreambuf_iterator<char>(is), {});
        it = input.begin();

        if (input.length() == 0) { throw NHXParserException("Error: empty input stream!\n"); }

        node_nothing(0, -1);
    }

    const AnnotatedTree& get_tree() const final { return tree; }
};
