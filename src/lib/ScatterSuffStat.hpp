#pragma once

#include "NodeMultivariateProcess.hpp"

class ScatterSuffStat : public SuffStat {
  public:
    //! Constructor (with only the tree given as argument)
    explicit ScatterSuffStat(const Tree &intree)
        : tree(intree), dimensions(1), scattermatrix(EMatrix(1, 1)){};

    ~ScatterSuffStat() = default;

    void Clear() {
        dimensions = 1;
        scattermatrix = EMatrix::Zero(1, 1);
    }

    void AddSuffStat(NodeMultivariateProcess const &nodeprocessses) {
        dimensions = nodeprocessses.GetDimensions();
        scattermatrix = EMatrix::Zero(dimensions, dimensions);

        for (Tree::BranchIndex branch = 0; branch < tree.nb_branches(); branch++) {
            EVector contrast = nodeprocessses.GetContrast(tree.node_index(branch));
            scattermatrix += contrast * contrast.transpose();
        }
    }

    //! return log p(Y | M) as a function of the covariance matrix M
    double GetLogProb(CovMatrix const &cov_matrix) const {
        EMatrix inv_cov = cov_matrix.inverse();
        double norm = cov_matrix.norm();
        inv_cov = scattermatrix * inv_cov;
        double trace = inv_cov.trace();
        return -(tree.nb_branches() * log(norm) + trace) / 2;
    }

    //! return log p(M | Y) propto p(Y | M) p(M) as a function of the covariance matrix M
    double GetLogPosterior(CovMatrix const &cov_matrix, int df, double kappa) const {
        EMatrix inv_cov = cov_matrix.inverse();
        double norm = cov_matrix.norm();

        EMatrix diag = EMatrix::Zero(dimensions, dimensions);
        for (int i{0}; i < dimensions; i++) { diag(i, i) = kappa; }

        inv_cov = (diag + scattermatrix) * inv_cov;
        double trace = inv_cov.trace();
        return -((df + dimensions + tree.nb_branches() + 1) * log(norm) + trace) / 2;
    }

  protected:
    const Tree &tree;
    int dimensions;
    EMatrix scattermatrix;
};