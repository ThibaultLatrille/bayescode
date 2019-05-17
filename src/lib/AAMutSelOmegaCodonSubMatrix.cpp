#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include <tuple>

void AAMutSelOmegaCodonSubMatrix::ComputeStationary() const {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < Nstate; i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i));
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    // double min = 1;
    for (int i = 0; i < Nstate; i++) { mStationary[i] /= total; }
}

void AAMutSelOmegaCodonSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (auto j : statespace->GetNeighbors(i)) {
        int pos = GetDifferingPosition(i, j);
        int a = GetCodonPosition(pos, i);
        int b = GetCodonPosition(pos, j);

        Q(i, j) = (*NucMatrix)(a, b);

        if (!Synonymous(i, j)) {
            double deltaS = GetLogFitness(GetCodonStateSpace()->Translation(j)) -
                            GetLogFitness(GetCodonStateSpace()->Translation(i));
            if ((fabs(deltaS)) < 1e-30) {
                Q(i, j) *= 1 + deltaS / 2;
            } else if (deltaS > 50) {
                Q(i, j) *= deltaS;
            } else if (deltaS < -50) {
                Q(i, j) = 0;
            } else {
                Q(i, j) *= deltaS / (1.0 - exp(-deltaS));
            }

            Q(i, j) *= GetOmega();
        }

        total += Q(i, j);

        assert(!std::isinf(Q(i, j)));
        assert(Q(i, j) >= 0);
    }

    Q(i, i) = -total;
    assert(total >= 0);
}

std::tuple<double, double> AAMutSelOmegaCodonSubMatrix::GetFlowDNDS() const {
    UpdateStationary();
    double totom = 0;
    double totweight = 0;
    for (int i = 0; i < Nstate; i++) {
        double weight = 0;
        double om = 0;
        for (auto j : statespace->GetNeighbors(i)) {
            if (!Synonymous(i, j)) {
                int pos = GetDifferingPosition(i, j);
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);

                double nucrate = (*NucMatrix)(a, b);

                double deltaS = GetLogFitness(GetCodonStateSpace()->Translation(j)) -
                                GetLogFitness(GetCodonStateSpace()->Translation(i));
                double pfix;
                if ((fabs(deltaS)) < 1e-30) {
                    pfix = 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    pfix = deltaS;
                } else if (deltaS < -50) {
                    pfix = 0;
                } else {
                    pfix = deltaS / (1.0 - exp(-deltaS));
                }

                om += nucrate * pfix;
                weight += nucrate;
            }
        }

        totom += mStationary[i] * om;
        totweight += mStationary[i] * weight;
    }
    return std::make_tuple(totom, totweight);
}

double AAMutSelOmegaCodonSubMatrix::GetPredictedDNDS() const {
    double dn = 0, dn0 = 0;
    std::tie(dn, dn0) = GetFlowDNDS();
    return dn / dn0;
}
