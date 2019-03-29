#include <cassert>

#include "AAMutSelNeCodonMatrixBidimArray.hpp"

MutSelNeCodonMatrixBidimArray::MutSelNeCodonMatrixBidimArray(
    const CodonStateSpace *incodonstatespace, const SubMatrix *innucmatrix,
    const Selector<std::vector<double>> *infitnessarray, const std::vector<double> &pop_size_array)
    : codonstatespace(incodonstatespace),
      nucmatrix(innucmatrix),
      aafitnessarray(infitnessarray),
      matrixbidimarray(pop_size_array.size(),
          std::vector<AAMutSelOmegaCodonSubMatrix *>(infitnessarray->GetSize(), nullptr)) {
    std::cout << GetNrow() << "\t" << GetNcol() << "\n";
    for (int i = 0; i < GetNrow(); i++) {
        for (int j = 0; j < GetNcol(); j++) {
            matrixbidimarray[i][j] = new AAMutSelOmegaCodonSubMatrix(
                codonstatespace, nucmatrix, infitnessarray->GetVal(j), 1.0, pop_size_array.at(i));
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

void MutSelNeCodonMatrixBidimArray::SetRowNe(int i, double Ne) {
    for (int j = 0; j < this->GetNcol(); j++) { (*this)(i, j).SetNe(Ne); }
}

void MutSelNeCodonMatrixBidimArray::SetNe(std::vector<double> const &Ne) {
    assert(GetNcol() > 0);
    assert(GetNrow() > 0);
    assert(GetNrow() == static_cast<int>(Ne.size()));
    for (int i = 0; i < this->GetNrow(); i++) { SetRowNe(i, Ne[i]); }
}

void MutSelNeCodonMatrixBidimArray::CorruptCodonMatrices() {
    for (int j = 0; j < GetNcol(); j++) { CorruptColCodonMatrices(j); }
}

void MutSelNeCodonMatrixBidimArray::CorruptColCodonMatrices(int j) {
    assert(GetNrow() > 0);
    for (int i = 0; i < GetNrow(); i++) { CorruptMatrix(i, j); }
}

void MutSelNeCodonMatrixBidimArray::CorruptRowCodonMatrices(int i) {
    assert(GetNcol() > 0);
    for (int j = 0; j < GetNcol(); j++) { CorruptMatrix(i, j); }
}

void MutSelNeCodonMatrixBidimArray::CorruptMatrix(int i, int j) { (*this)(i, j).CorruptMatrix(); }

void MutSelNeCodonMatrixBidimArray::UpdateCodonMatrices(const Selector<int> &occupancy) {
    for (int col = 0; col < GetNcol(); col++) {
        if (!occupancy.GetVal(col)) { (*this).CorruptColCodonMatrices(col); }
    }
}

AAMutSelNeCodonSubMatrixArray::AAMutSelNeCodonSubMatrixArray(
    const CodonStateSpace *incodonstatespace, const SubMatrix *innucmatrix,
    const Selector<std::vector<double>> *inaafitnessarray, double inne)
    : codonstatespace(incodonstatespace),
      nucmatrix(innucmatrix),
      aafitnessarray(inaafitnessarray),
      ne(inne),
      matrixarray(inaafitnessarray->GetSize()) {
    Create();
}

void AAMutSelNeCodonSubMatrixArray::UpdateCodonMatrices() {
    for (int i = 0; i < GetSize(); i++) { UpdateCodonMatrices(i); }
}

void AAMutSelNeCodonSubMatrixArray::UpdateCodonMatrices(int i) {
    (*this)[i].SetNe(ne);
    (*this)[i].CorruptMatrix();
}

void AAMutSelNeCodonSubMatrixArray::Create() {
    for (int i = 0; i < GetSize(); i++) {
        matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(
            codonstatespace, nucmatrix, aafitnessarray->GetVal(i), 1.0, ne);
    }
}

void AAMutSelNeCodonSubMatrixArray::Delete() {
    for (int i = 0; i < GetSize(); i++) { delete matrixarray[i]; }
}
