#include "ScaledMutationRateCompound.hpp"


NodeProcessScaledMutationRate::NodeProcessScaledMutationRate(double const &intheta_scale,
    NodeProcess const *innode_popsize, NodeProcess const *innode_rates, int Ntaxa)
    : theta_scale{intheta_scale}, node_popsize{innode_popsize}, node_rates{innode_rates} {
    theta.resize(Ntaxa);
}

double NodeProcessScaledMutationRate::GetTheta(int taxon) const {
    Tree::NodeIndex node = taxon_map->TaxonToNode(taxon);
    assert(node_rates->GetTree().is_leaf(node));
    return theta_scale * node_rates->GetExpVal(node) * node_popsize->GetExpVal(node);
}

void NodeProcessScaledMutationRate::Update() {
    for (size_t taxon = 0; taxon < theta.size(); taxon++) { theta[taxon] = GetTheta(taxon); }
}

BranchWiseProcessScaledMutationRate::BranchWiseProcessScaledMutationRate(
    double const &intheta_scale, LeafMultivariateProcess &inleaf_multivariate_process, int Ntaxa)
    : theta_scale{intheta_scale}, leaf_multivariate_process{inleaf_multivariate_process} {
    theta.resize(Ntaxa);
}

void BranchWiseProcessScaledMutationRate::SlidingTaxonMove(int taxon, int dimension, double m) {
    Tree::NodeIndex node = taxon_map->TaxonToNode(taxon);
    assert(leaf_multivariate_process.GetTree().is_leaf(node));
    leaf_multivariate_process[node](dimension) += m;
    theta[taxon] = theta_scale * leaf_multivariate_process.GetTheta(node);
}


void BranchWiseProcessScaledMutationRate::Update() {
    for (size_t taxon = 0; taxon < theta.size(); taxon++) {
        theta[taxon] =
            theta_scale * leaf_multivariate_process.GetTheta(taxon_map->TaxonToNode(taxon));
    }
}
