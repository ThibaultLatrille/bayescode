#include <cassert>

#include "AAMutSelNeCodonMatrixBidimArray.hpp"

MutSelNeCodonMatrixBidimArray::MutSelNeCodonMatrixBidimArray(const CodonStateSpace *codonstatespace,
    const SubMatrix *nucmatrix, const Selector<std::vector<double>> *fitnessarray,
    const std::vector<double> &pop_size_array, const std::vector<double> &omega_array)
    : matrixbidimarray(pop_size_array.size(),
          std::vector<AAMutSelOmegaCodonSubMatrix *>(fitnessarray->GetSize(), nullptr)) {
    std::cout << GetNrow() << "\t" << GetNcol() << "\n";
    for (int i = 0; i < GetNrow(); i++) {
        double omega = omega_array.empty() ? 1.0 : omega_array.at(i);
        for (int j = 0; j < GetNcol(); j++) {
            matrixbidimarray[i][j] = new AAMutSelOmegaCodonSubMatrix(
                codonstatespace, nucmatrix, fitnessarray->GetVal(j), omega, pop_size_array.at(i));
        }
    }
    assert(static_cast<int>(matrixbidimarray.size()) == GetNrow());
}

MutSelNeCodonMatrixBidimArray::~MutSelNeCodonMatrixBidimArray() {
    for (int i = 0; i < GetNrow(); i++) {
        for (int j = 0; j < GetNcol(); j++) { delete matrixbidimarray[i][j]; }
    }
}

const AAMutSelOmegaCodonSubMatrix &MutSelNeCodonMatrixBidimArray::GetVal(int i, int j) const {
    assert(0 <= i and i < GetNrow());
    assert(0 <= j and j < GetNcol());
    return *matrixbidimarray.at(i).at(j);
}

AAMutSelOmegaCodonSubMatrix &MutSelNeCodonMatrixBidimArray::operator()(int i, int j) {
    assert(0 <= i and i < GetNrow());
    assert(0 <= j and j < GetNcol());
    return *matrixbidimarray.at(i).at(j);
}

void MutSelNeCodonMatrixBidimArray::UpdateRowNe(int i, double Ne) {
    for (int j = 0; j < this->GetNcol(); j++) { (*this)(i, j).UpdateNe(Ne); }
}

void MutSelNeCodonMatrixBidimArray::UpdateRowOmega(int i, double omega) {
    for (int j = 0; j < this->GetNcol(); j++) { (*this)(i, j).UpdateOmega(omega); }
}

void MutSelNeCodonMatrixBidimArray::UpdateCodonMatricesNoFitnessRecomput() {
    for (int i = 0; i < this->GetNrow(); i++) {
        for (int j = 0; j < GetNcol(); j++) { (*this)(i, j).CorruptMatrixNoFitnessRecomput(); }
    }
}

void MutSelNeCodonMatrixBidimArray::UpdateColCodonMatrices(int j) {
    assert(GetNrow() > 0);
    for (int i = 0; i < GetNrow(); i++) { UpdateMatrix(i, j); }
}

void MutSelNeCodonMatrixBidimArray::UpdateMatrix(int i, int j) { (*this)(i, j).CorruptMatrix(); }

void MutSelNeCodonMatrixBidimArray::UpdateColCodonMatrices(const Selector<int> &occupancy) {
    for (int col = 0; col < GetNcol(); col++) {
        if (!occupancy.GetVal(col)) { (*this).UpdateColCodonMatrices(col); }
    }
}

AAMutSelNeCodonSubMatrixArray::AAMutSelNeCodonSubMatrixArray(const CodonStateSpace *codonstatespace,
    const SubMatrix *nucmatrix, const Selector<std::vector<double>> *aafitnessarray, double ne,
    double omega)
    : matrixarray(aafitnessarray->GetSize()) {
    for (int i = 0; i < aafitnessarray->GetSize(); i++) {
        matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(
            codonstatespace, nucmatrix, aafitnessarray->GetVal(i), omega, ne);
    }
}

void AAMutSelNeCodonSubMatrixArray::Delete() {
    for (int i = 0; i < GetSize(); i++) { delete matrixarray[i]; }
}
