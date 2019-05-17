#pragma once

#include "MultivariateProcess.hpp"
#include "SuffStat.hpp"

class ScatterSuffStat : public SuffStat {
  public:
    //! Constructor (with only the tree given as argument)
    explicit ScatterSuffStat(const Tree &intree);

    ~ScatterSuffStat() = default;

    void Clear();

    void AddSuffStat(NodeMultivariateProcess const &nodeprocessses);

    void SamplePrecisionMatrix(EMatrix &sampling_matrix, int df, double kappa) const;

  protected:
    const Tree &tree;
    int dimensions;
    EMatrix scattermatrix;
};