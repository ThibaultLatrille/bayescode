
#include "AAMutSelOmegaCodonSubMatrix.hpp"

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
    for (int j = 0; j < Nstate; j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);

                Q(i, j) = (*NucMatrix)(a, b);

                double deltaS = 0;
                if (!Synonymous(i, j)) {
                    deltaS = GetLogFitness(GetCodonStateSpace()->Translation(j)) -
                             GetLogFitness(GetCodonStateSpace()->Translation(i));
                }
                if ((fabs(deltaS)) < 1e-30) {
                    Q(i, j) *= 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    Q(i, j) *= deltaS;
                } else if (deltaS < -50) {
                    Q(i, j) = 0;
                } else {
                    Q(i, j) *= deltaS / (1.0 - exp(-deltaS));
                }
                if (!Synonymous(i, j)) { Q(i, j) *= GetOmega(); }
            } else {
                Q(i, j) = 0;
            }
            total += Q(i, j);

            if (std::isinf(Q(i, j))) {
                std::cerr << "Q matrix infinite: " << Q(i, j) << '\n';
                exit(1);
            }

            if (Q(i, j) < 0) {
                std::cerr << "Q matrix negative: " << Q(i, j) << '\n';
                exit(1);
            }
        }
    }

    Q(i, i) = -total;

    if (total < 0) {
        std::cerr << "negative rate away\n";
        exit(1);
    }
}

double AAMutSelOmegaCodonSubMatrix::GetPredictedDNDS() const {
    UpdateMatrix();
    double totom = 0;
    double totweight = 0;
    for (int i = 0; i < Nstate; i++) {
        double weight = 0;
        double om = 0;

        for (int j = 0; j < Nstate; j++) {
            if (i != j) {
                if (!Synonymous(i, j)) {
                    int pos = GetDifferingPosition(i, j);
                    if ((pos != -1) && (pos != 3)) {
                        int a = GetCodonPosition(pos, i);
                        int b = GetCodonPosition(pos, j);

                        double nucrate = (*NucMatrix)(a, b);

                        double deltaS = GetLogFitness(GetCodonStateSpace()->Translation(j)) -
                                        GetLogFitness(GetCodonStateSpace()->Translation(i));
                        double pfix = 1.0;
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
            }
        }

        totom += mStationary[i] * om;
        totweight += mStationary[i] * weight;
    }
    return totom / totweight;
}
