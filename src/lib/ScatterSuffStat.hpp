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
    double GetLogProb(PrecisionMatrix const &precision_matrix) const {
        double det = precision_matrix.determinant();
        EMatrix mul = scattermatrix * precision_matrix;
        double trace = mul.trace();
        return 0.5 * tree.nb_branches() * log(det) - 0.5 * trace;
    }

    //! return log p(M | Y) propto p(Y | M) p(M) as a function of the covariance matrix M
    double GetLogPosterior(PrecisionMatrix const &precision_matrix, int df, double kappa) const {
        double det = precision_matrix.determinant();

        EMatrix diag = EMatrix::Zero(dimensions, dimensions);
        for (int i{0}; i < dimensions; i++) { diag(i, i) = kappa; }

        EMatrix mul = (diag + scattermatrix) * precision_matrix;
        double trace = mul.trace();
        return 0.5 * (df + dimensions + tree.nb_branches() + 1) * log(det) - 0.5 * trace;
    }

  protected:
    const Tree &tree;
    int dimensions;
    EMatrix scattermatrix;
};