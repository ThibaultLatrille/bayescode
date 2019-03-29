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

double NodeProcess::GetExpVal(Tree::NodeIndex node) const { return exp(GetVal(node)); }

void NodeProcess::SlidingMove(Tree::NodeIndex node, double m) {
    node_multivariate[node](dimension) += m;
}

void NodeProcess::SlidingMove(double m) {
    for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(GetTree().nb_nodes()); node++) {
        SlidingMove(node, m);
    }
}

BranchProcess::BranchProcess(const NodeProcess &innodeprocess)
    : SimpleBranchArray<double>(innodeprocess.GetTree()), nodeprocess{innodeprocess} {
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
    (*this)[this->GetTree().branch_index(node)] =
        (nodeprocess.GetExpVal(parent) + nodeprocess.GetExpVal(node)) / 2;
    // For geodesic, use :
    // = (nodeprocess.GetExpVal(parent) - nodeprocess.GetExpVal(node)) /
    // (nodeprocess->GetVal(parent) - nodeprocess->GetVal(node));
}

BranchWiseMultivariateProcess::BranchWiseMultivariateProcess(
    const Chronogram &inchrono, const EMatrix &inprecision_matrix, int indimensions)
    : SimpleBranchArray<EVector>(inchrono.GetTree()),
      chronogram(inchrono),
      dimensions(indimensions),
      precision_matrix(inprecision_matrix) {
    root_process = EVector::Zero(dimensions);
    for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
        (*this)[branch] = EVector::Zero(dimensions);
    }
}

double BranchWiseMultivariateProcess::GetBranchLogProb(Tree::BranchIndex branch) const {
    return Random::logNormalDensity(GetContrast(branch), precision_matrix);
}

double BranchWiseMultivariateProcess::GetLocalBranchLogProb(Tree::BranchIndex branch) const {
    double tot = GetBranchLogProb(branch);
    for (Tree::NodeIndex const &child : GetTree().children(GetTree().node_index(branch))) {
        tot += GetBranchLogProb(GetTree().branch_index(child));
    }
    return tot;
}

EVector BranchWiseMultivariateProcess::GetContrast(Tree::BranchIndex branch) const {
    double distance = chronogram.GetVal(branch) / 2;
    Tree::NodeIndex node = GetTree().node_index(branch);
    assert(!GetTree().is_root(node));
    Tree::NodeIndex parent_node = GetTree().parent(node);
    if (GetTree().is_root(parent_node)) {
        return (root_process - this->GetVal(branch)) / sqrt(distance);
    } else {
        Tree::BranchIndex parent_branch = GetTree().branch_index(parent_node);
        distance += chronogram.GetVal(parent_branch) / 2;
        return (this->GetVal(parent_branch) - this->GetVal(branch)) / sqrt(distance);
    }
}

double BranchWiseMultivariateProcess::GetLocalNodeLogProb(Tree::NodeIndex node) const {
    double tot = 0;
    if (!GetTree().is_root(node)) { tot += GetBranchLogProb(GetTree().branch_index(node)); }
    for (Tree::NodeIndex const &child : GetTree().children(node)) {
        tot += GetBranchLogProb(GetTree().branch_index(child));
    }
    return tot;
}

double BranchWiseMultivariateProcess::GetLogProb() const {
    double tot = 0;
    for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
        tot += GetBranchLogProb(branch);
    }
    return tot;
}

LeafMultivariateProcess::LeafMultivariateProcess(
    const BranchWiseMultivariateProcess &inbranchwiseprocess, TaxonSet const &taxon)
    : SimpleNodeArray<EVector>(inbranchwiseprocess.GetTree()),
      branchwiseprocess(inbranchwiseprocess) {
    index_table = taxon.get_index_table(&branchwiseprocess.GetTree());
    for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
        (*this)[branch] = EVector::Zero(branchwiseprocess.GetDimensions());
    }
}

double LeafMultivariateProcess::GetTaxonLogProb(int taxon) const {
    return GetLogProb(index_table[taxon]);
}

double LeafMultivariateProcess::GetLogProb(Tree::NodeIndex node) const {
    assert(GetTree().is_leaf(node));
    return Random::logNormalDensity(GetContrast(node), branchwiseprocess.precision_matrix);
}

EVector LeafMultivariateProcess::GetContrast(Tree::NodeIndex node) const {
    Tree::BranchIndex parent_branch = GetTree().branch_index(node);
    double distance = branchwiseprocess.chronogram.GetVal(parent_branch) / 2;
    return (branchwiseprocess.GetVal(parent_branch) - this->GetVal(node)) / sqrt(distance);
}

double LeafMultivariateProcess::GetTheta(Tree::NodeIndex node) const {
    assert(GetTree().is_leaf(node));
    return exp(this->GetVal(node).sum());
}

double LeafMultivariateProcess::GetLogProb() const {
    double tot = 0;
    for (Tree::NodeIndex node = 0; node < static_cast<Tree::NodeIndex>(GetTree().nb_nodes());
         node++) {
        if (tree.is_leaf(node)) { tot += GetLogProb(node); }
    }
    return tot;
}

BranchWiseProcess::BranchWiseProcess(
    BranchWiseMultivariateProcess &inbranchwise_multivariate, int indimension)
    : SimpleBranchArray<double>(inbranchwise_multivariate.GetTree()),
      dimension(indimension),
      branchwise_multivariate(inbranchwise_multivariate) {
    Update();
}

void BranchWiseProcess::SlidingMove(Tree::BranchIndex branch, double m) {
    branchwise_multivariate[branch](dimension) += m;
    this->array[branch] = exp(branchwise_multivariate[branch](dimension));
}

void BranchWiseProcess::SlidingRootMove(double m) {
    branchwise_multivariate.root_process(dimension) += m;
    root_value = exp(branchwise_multivariate.root_process(dimension));
}

void BranchWiseProcess::Update() {
    root_value = exp(branchwise_multivariate.root_process(dimension));
    for (Tree::BranchIndex branch = 0; branch < GetTree().nb_branches(); branch++) {
        this->array[branch] = exp(branchwise_multivariate[branch](dimension));
    }
}
