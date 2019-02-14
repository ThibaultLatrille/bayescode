#pragma once

#include "NodeMultivariateProcess.hpp"

/**
 * \brief A ScaledMutationRate that returns the scaled mutation rate (theta=4*Ne*u) as a compound
 * parameter given Ne (population size) and u (mutation rate per generation) of a given taxon.
 *
 * Useful for implementing models assuming each taxon has its own Ne (population size) and u
 * (mutation rate per generation). Moreover Ne and u are assumed to follow a brownian process
 * (NodeProcess).
 */
class CompoundScaledMutationRate : public ScaledMutationRate {
  public:
    //! \brief Constructor, taking as arguments the Nbr of taxa, the scaling factor, the mutation
    //! rate (a NodeProcess) and the population size (a NodeProcess).
    CompoundScaledMutationRate(int Ntaxa, double const &intheta_scale,
        NodeProcess const *innode_rates, NodeProcess const *innode_popsize, TaxonSet const &taxon)
        : theta_scale{intheta_scale}, node_rates{innode_rates}, node_popsize{innode_popsize} {
        reverse_index_table = taxon.get_reverse_index_table(&node_rates->GetTree());
    };

    ~CompoundScaledMutationRate() override = default;

    double GetTheta(int taxon) const override {
        Tree::NodeIndex node = reverse_index_table[taxon];
        return theta_scale * node_rates->GetExpVal(node) * node_popsize->GetExpVal(node);
    }

  private:
    double const &theta_scale;
    NodeProcess const *node_rates;
    NodeProcess const *node_popsize;

    std::vector<Tree::NodeIndex> reverse_index_table;
};
