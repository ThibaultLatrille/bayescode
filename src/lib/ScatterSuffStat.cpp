#include "ScatterSuffStat.hpp"

ScatterSuffStat::ScatterSuffStat(const Tree &intree)
    : tree(intree), dimensions(1), scattermatrix(EMatrix(1, 1)) {}

void ScatterSuffStat::Clear() {
    dimensions = 1;
    scattermatrix = EMatrix::Zero(1, 1);
}

void ScatterSuffStat::AddSuffStat(NodeMultivariateProcess const &nodeprocessses) {
    dimensions = nodeprocessses.GetDimensions();
    scattermatrix = EMatrix::Zero(dimensions, dimensions);

    for (Tree::BranchIndex branch = 0; branch < tree.nb_branches(); branch++) {
        EVector contrast = nodeprocessses.GetContrast(tree.node_index(branch));
        scattermatrix += contrast * contrast.transpose();
    }
}

void ScatterSuffStat::SamplePrecisionMatrix(
    PrecisionMatrix &precision_matrix, PriorCovariance const &prior_cov_matrix) const {
    precision_matrix.setZero();

    EMatrix cov = (prior_cov_matrix.GetPriorCovarianceMatrix() + scattermatrix).inverse();
    Eigen::SelfAdjointEigenSolver<EMatrix> e_solver(cov);
    EMatrix A = e_solver.eigenvectors() * e_solver.eigenvalues().cwiseSqrt().asDiagonal();
    assert((A * A.transpose() - cov).norm() < 1e-8);
    EVector sampled_vector = EVector::Zero(dimensions);

    int nbr_samples = prior_cov_matrix.GetDoF() + tree.nb_branches();
    for (int i = 0; i < nbr_samples; i++) {
        for (int dim = 0; dim < dimensions; dim++) { sampled_vector(dim) = Random::sNormal(); }
        sampled_vector = A * sampled_vector;
        precision_matrix += sampled_vector * sampled_vector.transpose();
    }
}
