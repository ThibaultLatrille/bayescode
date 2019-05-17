#pragma once

#include "MultivariateProcess.hpp"
#include "ScaledMutationRate.hpp"
#include "TaxonMapping.hpp"

/**
 * \brief A ScaledMutationRate that returns the scaled mutation rate (theta=4*Ne*u) as a compound
 * parameter given Ne (population size) and u (mutation rate per generation) of a given taxon.
 *
 * Useful for implementing models assuming each taxon has its own Ne (population size) and u
 * (mutation rate per generation). Moreover Ne and u are assumed to follow a brownian process
 * (NodeProcess).
 */
class NodeProcessScaledMutationRate : public ScaledMutationRate {
  public:
    //! \brief Constructor, taking as arguments the Nbr of taxa, the scaling factor, the mutation
    //! rate (a NodeProcess) and the population size (a NodeProcess).
    NodeProcessScaledMutationRate(double const &intheta_scale, NodeProcess const *innode_popsize,
        NodeProcess const *innode_rates, int Ntaxa);

    ~NodeProcessScaledMutationRate() override = default;

    void SetTaxonMap(TaxonMap const *intaxon_map) { taxon_map = intaxon_map; }

    double GetTheta(int taxon) const override;

    double GetNodeTheta(Tree::NodeIndex node) const;

    double &operator[](int taxon) { return theta[taxon]; }

    void Update();

  private:
    double const &theta_scale;
    NodeProcess const *node_popsize;
    NodeProcess const *node_rates;
    std::vector<double> theta;

    TaxonMap const *taxon_map{nullptr};
};