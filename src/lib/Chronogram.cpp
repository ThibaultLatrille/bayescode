#include <cassert>
#include "Chronogram.hpp"

NodeAges::NodeAges(const Tree &intree) : SimpleNodeArray<double>(intree) { UniformSample(); }

void NodeAges::UniformSample() {
    double eccentricity = EccentricityRecursive(GetTree().root());

    for (Tree::NodeIndex node{0}; node < GetNnode(); node++) {
        // renormalize ages so that root age is 1 by definition
        (*this)[node] /= eccentricity;
    }
    Check();
}

double NodeAges::EccentricityRecursive(Tree::NodeIndex node) {
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
}

void NodeAges::Check() const {
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

void NodeAges::SlidingMove(Tree::NodeIndex node, double scale) {
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

Chronogram::Chronogram(const NodeAges &innodeages)
        : SimpleBranchArray<double>(innodeages.GetTree()), nodeages(innodeages) {
    Update();
}

void Chronogram::Update() {
    for (Tree::NodeIndex node{0}; node < Tree::NodeIndex(this->GetTree().nb_nodes()); node++) {
        if (!this->GetTree().is_root(node)) {
            this->UpdateBranch(this->GetTree().parent(node), node);
        }
    }
}

void Chronogram::UpdateLocal(Tree::NodeIndex node) {
    // update all branch lengths around this node

    // for all children
    for (auto const& child : this->GetTree().children(node)) { UpdateBranch(node, child); }

    // for the branch attached to the node
    if (!this->GetTree().is_root(node)) { UpdateBranch(this->GetTree().parent(node), node); }
}

void Chronogram::UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node) {
    (*this)[this->GetTree().branch_index(node)] =
            nodeages.GetVal(parent) - nodeages.GetVal(node);
}
