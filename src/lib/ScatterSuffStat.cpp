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

void ScatterSuffStat::SamplePrecisionMatrix(EMatrix &sampling_matrix, int df, double kappa) const {
    sampling_matrix.setZero();

    EMatrix prior_matrix = EMatrix::Identity(dimensions, dimensions) * kappa;
    EMatrix precision_matrix = (prior_matrix + scattermatrix).inverse();

    Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(precision_matrix);
    EMatrix transform =
        eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();
    EVector sampled_vector = EVector::Zero(dimensions);

    int nbr_samples = df + tree.nb_branches();
    for (int i = 0; i < nbr_samples; i++) {
        for (int dim = 0; dim < dimensions; dim++) { sampled_vector(dim) = Random::sNormal(); }
        sampled_vector = transform * sampled_vector;
        sampling_matrix += sampled_vector * sampled_vector.transpose();
    }
}
