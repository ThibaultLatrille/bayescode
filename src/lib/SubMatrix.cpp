#include "SubMatrix.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

int SubMatrix::nuni = 0;
int SubMatrix::nunimax = 0;
int SubMatrix::nunisubcount = 0;
int SubMatrix::diagcount = 0;
double SubMatrix::diagerr = 0;

double SubMatrix::nz = 0;
double SubMatrix::meanz = 0;
double SubMatrix::maxz = 0;

const int witheigen = 1;

// ---------------------------------------------------------------------------
//     SubMatrix()
// ---------------------------------------------------------------------------
SubMatrix::SubMatrix(int inNstate, bool innormalise) : Nstate(inNstate), normalise(innormalise) {
    ndiagfailed = 0;
    Create();
}

void SubMatrix::Create() {
    Q = EMatrix::Zero(Nstate, Nstate);
    u = EMatrix(Nstate, Nstate);
    invu = EMatrix(Nstate, Nstate);
    v = EVector(Nstate);
    vi = EVector(Nstate);
    mStationary = EVector(Nstate);

    ptrQ = nullptr;
    ptru = nullptr;
    ptrinvu = nullptr;
    ptrv = nullptr;
    ptrStationary = nullptr;

    if (!witheigen) {
        ptrQ = new double *[Nstate];
        for (int i = 0; i < Nstate; i++) { ptrQ[i] = new double[Nstate]; }

        ptru = new double *[Nstate];
        for (int i = 0; i < Nstate; i++) { ptru[i] = new double[Nstate]; }

        ptrinvu = new double *[Nstate];
        for (int i = 0; i < Nstate; i++) { ptrinvu[i] = new double[Nstate]; }

        ptrv = new double[Nstate];
        ptrStationary = new double[Nstate];
    }

    UniMu = 1;
    mPow = new double **[UniSubNmax];
    for (int n = 0; n < UniSubNmax; n++) { mPow[n] = nullptr; }

    flagarray = new bool[Nstate];
    diagflag = false;
    statflag = false;
    for (int i = 0; i < Nstate; i++) { flagarray[i] = false; }
    powflag = false;
}

// ---------------------------------------------------------------------------
//     ~SubMatrix()
// ---------------------------------------------------------------------------

SubMatrix::~SubMatrix() {
    if (mPow != nullptr) {
        for (int n = 0; n < UniSubNmax; n++) {
            if (mPow[n] != nullptr) {
                for (int i = 0; i < Nstate; i++) { delete[] mPow[n][i]; }
                delete[] mPow[n];
            }
        }
        delete[] mPow;
    }
    delete[] flagarray;
}

// ---------------------------------------------------------------------------
//     void ScalarMul()
// ---------------------------------------------------------------------------

void SubMatrix::ScalarMul(double e) {
    for (int i = 0; i < Nstate; i++) {
        v[i] *= e;
        vi[i] *= e;
    }
    UniMu *= e;
}

// ---------------------------------------------------------------------------
//     Diagonalise()
// ---------------------------------------------------------------------------

int SubMatrix::Diagonalise() const {
    if (witheigen) {
        EigenDiagonalise();
    } else {
        OldDiagonalise();
    }
    return 0;
}

int SubMatrix::EigenDiagonalise() const {
    if (!ArrayUpdated()) { UpdateMatrix(); }

    diagcount++;
    auto &stat = GetStationary();

    EMatrix a(Nstate, Nstate);

    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) { a(i, j) = Q(i, j) * sqrt(stat[i] / stat[j]); }
    }

    solver.compute(a);
    v = solver.eigenvalues().real();
    vi = solver.eigenvalues().imag();
    u = solver.eigenvectors().real();

    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) { invu(i, j) = u(j, i) * sqrt(stat[j]); }
    }
    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) { u(i, j) /= sqrt(stat[i]); }
    }

    diagflag = true;
    double err = CheckDiag();
    if (diagerr < err) { diagerr = err; }
    return 0;
}

double SubMatrix::CheckDiag() const {
    EMatrix tmp(Nstate, Nstate);
    EMatrix D(Nstate, Nstate);
    EMatrix Q2(Nstate, Nstate);

    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) { D(i, j) = 0; }
    }
    for (int i = 0; i < Nstate; i++) { D(i, i) = v[i]; }

    tmp = D * invu;
    Q2 = u * tmp;

    double max = 0;
    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) {
            double temp = fabs(Q2(i, j) - Q(i, j));
            if (max < temp) { max = temp; }
        }
    }
    return max;
}

// ---------------------------------------------------------------------------
//     ComputeRate()
// ---------------------------------------------------------------------------

double SubMatrix::GetRate() const {
    if (!ArrayUpdated()) {
        UpdateStationary();
        for (int k = 0; k < Nstate; k++) { ComputeArray(k); }
    }
    for (int k = 0; k < Nstate; k++) { flagarray[k] = true; }
    double norm = 0;
    for (int i = 0; i < Nstate - 1; i++) {
        for (int j = i + 1; j < Nstate; j++) { norm += mStationary[i] * Q(i, j); }
    }
    return 2 * norm;
}

// ---------------------------------------------------------------------------
//     Update()
// ---------------------------------------------------------------------------

void SubMatrix::UpdateMatrix() const {
    UpdateStationary();
    for (int k = 0; k < Nstate; k++) { ComputeArray(k); }
    for (int k = 0; k < Nstate; k++) { flagarray[k] = true; }
    if (isNormalised()) { Normalise(); }
}

// ---------------------------------------------------------------------------
//     Normalise()
// ---------------------------------------------------------------------------

void SubMatrix::Normalise() const {
    double norm = GetRate();
    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) { Q(i, j) /= norm; }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     Powers
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void SubMatrix::ActivatePowers() const {
    if (!powflag) {
        if (!ArrayUpdated()) { UpdateMatrix(); }

        UniMu = 0;
        for (int i = 0; i < Nstate; i++) {
            if (UniMu < fabs(Q(i, i))) { UniMu = fabs(Q(i, i)); }
        }

        CreatePowers(0);
        for (int i = 0; i < Nstate; i++) {
            for (int j = 0; j < Nstate; j++) { mPow[0][i][j] = 0; }
        }
        for (int i = 0; i < Nstate; i++) { mPow[0][i][i] = 1; }
        for (int i = 0; i < Nstate; i++) {
            for (int j = 0; j < Nstate; j++) {
                mPow[0][i][j] += Q(i, j) / UniMu;
                if (mPow[0][i][j] < 0) {
                    cerr << "error in SubMatrix::ComputePowers: negative prob : ";
                    cerr << i << '\t' << j << '\t' << mPow[0][i][j] << '\n';
                    cerr << "Nstate : " << Nstate << '\n';
                    exit(1);
                }
            }
        }
        npow = 1;
        powflag = true;
    }
}

void SubMatrix::InactivatePowers() const {
    if (powflag) {
        for (int n = 0; n < UniSubNmax; n++) {
            if (mPow[n] != nullptr) {
                for (int i = 0; i < Nstate; i++) { delete[] mPow[n][i]; }
                delete[] mPow[n];
                mPow[n] = nullptr;
            }
        }
        nunimax += npow;
        nuni++;

        npow = 0;
        powflag = false;
    }
}

void SubMatrix::CreatePowers(int n) const {
    if (mPow[n] == nullptr) {
        mPow[n] = new double *[Nstate];
        for (int i = 0; i < Nstate; i++) { mPow[n][i] = new double[Nstate]; }
    }
}

double SubMatrix::GetUniformizationMu() const {
    if (!powflag) { ActivatePowers(); }
    return UniMu;
}

double SubMatrix::Power(int n, int i, int j) const {
    if (!powflag) { ActivatePowers(); }
    if (n == 0) { return static_cast<double>(i == j); }
    if (n > UniSubNmax) { return Stationary(j); }
    if (n > npow) { ComputePowers(n); }
    return mPow[n - 1][i][j];
}

void SubMatrix::ComputePowers(int N) const {
    if (!powflag) { ActivatePowers(); }
    if (N > npow) {
        for (int n = npow; n < N; n++) {
            CreatePowers(n);
            for (int i = 0; i < Nstate; i++) {
                for (int j = 0; j < Nstate; j++) {
                    double &t = mPow[n][i][j];
                    t = 0;
                    for (int k = 0; k < Nstate; k++) { t += mPow[n - 1][i][k] * mPow[0][k][j]; }
                }
            }
        }
        npow = N;
    }
}

void SubMatrix::ToStream(ostream &os) const {
    os << GetNstate() << '\n';
    os << "stationaries: \n";
    for (int i = 0; i < GetNstate(); i++) { os << Stationary(i) << '\t'; }
    os << '\n';

    os << "rate matrix\n";
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            os << Q(i, j) << '\t';
            // os << Q[i][j] << '\t';
        }
        os << '\n';
    }

    os << '\n';
    for (int i = 0; i < GetNstate(); i++) { os << v[i] << '\t'; }
    os << '\n';
}

void SubMatrix::CheckReversibility() const {
    double max = 0;
    int imax = 0;
    int jmax = 0;
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = i + 1; j < GetNstate(); j++) {
            double tmp = fabs(Stationary(i) * Q(i, j) - Stationary(j) * Q(j, i));
            if (max < tmp) {
                max = tmp;
                imax = i;
                jmax = j;
            }
        }
    }
    if (max > 1e-6) {
        cerr << "max irreversibility: " << max << '\n';
        cerr << imax << '\t' << jmax << '\t' << Stationary(imax) << '\t' << Q(imax, jmax) << '\t'
             << Stationary(jmax) << '\t' << Q(jmax, imax) << '\n';
        exit(1);
    }
    cerr << "max rev: " << max << '\n';
}
