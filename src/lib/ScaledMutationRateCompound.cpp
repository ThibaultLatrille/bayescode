#include "ScaledMutationRateCompound.hpp"


NodeProcessScaledMutationRate::NodeProcessScaledMutationRate(double const &intheta_scale,
    NodeProcess const *innode_popsize, NodeProcess const *innode_rates, int Ntaxa)
    : theta_scale{intheta_scale}, node_popsize{innode_popsize}, node_rates{innode_rates} {
    theta.resize(Ntaxa);
}

double NodeProcessScaledMutationRate::GetTheta(int taxon) const {
    return GetNodeTheta(taxon_map->TaxonToNode(taxon));
}

double NodeProcessScaledMutationRate::GetNodeTheta(Tree::NodeIndex node) const {
    assert(node_rates->GetTree().is_leaf(node));
    return theta_scale * exp(node_rates->GetVal(node) + node_popsize->GetVal(node));
}

void NodeProcessScaledMutationRate::Update() {
    for (size_t taxon = 0; taxon < theta.size(); taxon++) { theta[taxon] = GetTheta(taxon); }
}