#include "CodonSubMatrix.hpp"
using namespace std;

void MGCodonSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (auto j : statespace->GetNeighbors(i)) {
        int pos = GetDifferingPosition(i, j);
        int a = GetCodonPosition(pos, i);
        int b = GetCodonPosition(pos, j);

        assert(a != b);
        Q(i, j) = (*NucMatrix)(a, b);
        total += Q(i, j);
    }
    Q(i, i) = -total;
    assert(total >= 0);
}

void MGCodonSubMatrix::ComputeStationary() const {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i));
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) { mStationary[i] /= total; }
}

void MGOmegaCodonSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (auto j : statespace->GetNeighbors(i)) {
        int pos = GetDifferingPosition(i, j);
        int a = GetCodonPosition(pos, i);
        int b = GetCodonPosition(pos, j);

        assert(a != b);
        Q(i, j) = (*NucMatrix)(a, b);
        if (!Synonymous(i, j)) { Q(i, j) *= GetOmega(); }

        total += Q(i, j);
    }
    Q(i, i) = -total;
    assert(total >= 0);
}