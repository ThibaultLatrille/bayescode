#include <cassert>
#include <cmath>

#include "MultivariateProcess.hpp"

NodeMultivariateProcess::NodeMultivariateProcess(
    const Chronogram &inchrono, const EMatrix &inprecision_matrix, int indimensions)
    : SimpleNodeArray<EVector>(inchrono.GetTree()),
      chronogram(inchrono),
      dimensions(indimensions),
      precision_matrix(inprecision_matrix) {
    for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
        (*this)[node] = EVector::Zero(dimensions);
    }
}

double NodeMultivariateProcess::GetSigma(int dimension) const {
    return sqrt(1.0 / precision_matrix(dimension, dimension));
}

double NodeMultivariateProcess::GetLogProb(Tree::NodeIndex node) const {
    if (GetTree().is_root(node)) {
        return 0.0;
    } else {
        return Random::logNormalDensity(GetContrast(node), precision_matrix);
    }
}

EVector NodeMultivariateProcess::GetContrast(Tree::NodeIndex node) const {
    assert(!GetTree().is_root(node));
    return (this->GetVal(node) - this->GetVal(GetTree().parent(node))) /
           sqrt(chronogram.GetVal(GetTree().branch_index(node)));
}

double NodeMultivariateProcess::GetLocalLogProb(Tree::NodeIndex node) const {
    double tot = GetLogProb(node);
    for (auto const &child : GetTree().children(node)) { tot += GetLogProb(child); }
    return tot;
}

double NodeMultivariateProcess::GetLogProb() const {
    double tot = 0;
    for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
        if (!GetTree().is_root(node)) { tot += GetLogProb(node); }
    }
    return tot;
}

NodeProcess::NodeProcess(NodeMultivariateProcess &innode_multivariate, int indimension)
    : dimension(indimension), node_multivariate(innode_multivariate) {}

double NodeProcess::GetSigma() const { return node_multivariate.GetSigma(dimension); }

const Tree &NodeProcess::GetTree() const { return node_multivariate.GetTree(); }

double NodeProcess::GetVal(Tree::NodeIndex node) const {
    return node_multivariate.GetVal(node)(dimension);
}

double &NodeProcess::operator[](Tree::NodeIndex node) { return node_multivariate[node](dimension); }

void NodeProcess::SlidingMove(Tree::NodeIndex node, double m) {
    node_multivariate[node](dimension) += m;
}

void NodeProcess::SlidingMove(double m) {
    for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
        SlidingMove(node, m);
    }
}

BranchProcess::BranchProcess(const NodeProcess &innodeprocess, bool arithmetic)
    : SimpleBranchArray<double>(innodeprocess.GetTree()),
      nodeprocess{innodeprocess},
      arithmetic{arithmetic} {
    Update();
}

void BranchProcess::Update() {
    for (Tree::NodeIndex node{0}; node < Tree::NodeIndex(this->GetTree().nb_nodes()); node++) {
        if (!this->GetTree().is_root(node)) {
            this->UpdateBranch(this->GetTree().parent(node), node);
        }
    }
}

void BranchProcess::UpdateLocal(Tree::NodeIndex node) {
    // update all branch lengths around this node

    // for all children
    for (auto const &child : this->GetTree().children(node)) { UpdateBranch(node, child); }

    // for the branch attached to the node
    if (!this->GetTree().is_root(node)) { UpdateBranch(this->GetTree().parent(node), node); }
}

void BranchProcess::UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node) {
    double xup = nodeprocess.GetVal(parent);
    double x = nodeprocess.GetVal(node);
    if (arithmetic) {
        (*this)[this->GetTree().branch_index(node)] = (exp(xup) + exp(x)) / 2;
    } else {
        if (abs(x - xup) < 1e-12) {
            (*this)[this->GetTree().branch_index(node)] = exp(x);
        } else {
            (*this)[this->GetTree().branch_index(node)] = (exp(xup) - exp(x)) / (xup - x);
        }
        assert((*this)[this->GetTree().branch_index(node)] >= 0.0);
    }
}