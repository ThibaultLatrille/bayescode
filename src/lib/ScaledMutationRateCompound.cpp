#include "ScaledMutationRateCompound.hpp"

NodeProcessScaledMutationRate::NodeProcessScaledMutationRate(int Ntaxa, double const &intheta_scale,
    NodeProcess const *innode_rates, NodeProcess const *innode_popsize, TaxonSet const &taxon)
    : theta_scale{intheta_scale}, node_rates{innode_rates}, node_popsize{innode_popsize} {
    reverse_index_table = taxon.get_reverse_index_table(&node_rates->GetTree());
    theta.resize(taxon.GetNtaxa());
}

double NodeProcessScaledMutationRate::GetTheta(int taxon) const {
    Tree::NodeIndex node = reverse_index_table[taxon];
    assert(node_rates->GetTree().is_leaf(node));
    return theta_scale * node_rates->GetExpVal(node) * node_popsize->GetExpVal(node);
}

void NodeProcessScaledMutationRate::Update() {
    for (size_t taxon = 0; taxon < theta.size(); taxon++) { theta[taxon] = GetTheta(taxon); }
}

BranchWiseProcessScaledMutationRate::BranchWiseProcessScaledMutationRate(int Ntaxa,
    double const &intheta_scale, LeafMultivariateProcess &inleaf_multivariate_process,
    TaxonSet const &taxon)
    : theta_scale{intheta_scale}, leaf_multivariate_process{inleaf_multivariate_process} {
    reverse_index_table = taxon.get_reverse_index_table(&leaf_multivariate_process.GetTree());
    theta.resize(taxon.GetNtaxa());
}

void BranchWiseProcessScaledMutationRate::SlidingTaxonMove(int taxon, int dimension, double m) {
    Tree::NodeIndex node = reverse_index_table[taxon];
    assert(leaf_multivariate_process.GetTree().is_leaf(node));
    leaf_multivariate_process[node](dimension) += m;
    theta[taxon] = theta_scale * leaf_multivariate_process.GetTheta(node);
}


void BranchWiseProcessScaledMutationRate::Update() {
    for (size_t taxon = 0; taxon < theta.size(); taxon++) {
        theta[taxon] = theta_scale * leaf_multivariate_process.GetTheta(reverse_index_table[taxon]);
    }
}
