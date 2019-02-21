#pragma once

#include <cassert>
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "Array.hpp"
#include "BidimArray.hpp"
#include "IIDDirichlet.hpp"

/**
 * \brief A BidimArray of AAMutSelOmegaCodonSubMatrix
 *
 * The constructor takes as arguments a single nucleotide matrix, an array of amino-acid
 * fitness profiles (across sites) and an array of population size (across branches). It then
 * constructs a 2-dimensional-array of AAMutSelOmegaCodonSubMatrix where the number of columns
 * is the size of the fitness profiles array (number of sites or component) and the number of rows
 * equal the size of the population size array (number of branches or condition).
 */

class MutSelNeCodonMatrixBidimArray : public BidimArray<SubMatrix>,
                                      public BidimArray<AAMutSelOmegaCodonSubMatrix> {
  public:
    //! constructor parameterized by a codon state space, a single nucleotide matrix, an array of
    //! fitness profiles and an array of population size.
    MutSelNeCodonMatrixBidimArray(const CodonStateSpace *incodonstatespace,
        const SubMatrix *innucmatrix, const Selector<std::vector<double>> *infitnessarray,
        const std::vector<double> *inNe)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          aafitnessarray(infitnessarray),
          popsizearray(inNe),
          matrixbidimarray(inNe->size(),
              std::vector<AAMutSelOmegaCodonSubMatrix *>(infitnessarray->GetSize(), nullptr)) {
        std::cout << GetNrow() << "\t" << GetNcol() << "\n";
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                matrixbidimarray[i][j] = new AAMutSelOmegaCodonSubMatrix(codonstatespace, nucmatrix,
                    infitnessarray->GetVal(j), 1.0, popsizearray->at(i));
            }
        }
        assert(static_cast<int>(matrixbidimarray.size()) == GetNrow());
    }

    ~MutSelNeCodonMatrixBidimArray() override {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { delete matrixbidimarray[i][j]; }
        }
    }

    //! return the number of rows (number of branches)
    int GetNrow() const override { return popsizearray->size(); }
    //! return the number of columns (number of sites)
    int GetNcol() const override { return aafitnessarray->GetSize(); }

    //! const access to the matrix for row (branch) i and column (site) j
    const AAMutSelOmegaCodonSubMatrix &GetVal(int i, int j) const override {
        assert(0 <= i and i < GetNrow());
        assert(0 <= j and j < GetNcol());
        return *matrixbidimarray.at(i).at(j);
    }

    //! non const access to the matrix for row (branch) i and column (site) j
    AAMutSelOmegaCodonSubMatrix &operator()(int i, int j) override {
        assert(0 <= i and i < GetNrow());
        assert(0 <= j and j < GetNcol());
        return *matrixbidimarray.at(i).at(j);
    }

    //! Set ne for row (branch) i
    void SetRowNe(int i, double Ne) {
        for (int j = 0; j < this->GetNcol(); j++) { (*this)(i, j).SetNe(Ne); }
    }

    void SetNe(std::vector<double> &Ne) {
        assert(GetNcol() > 0);
        assert(GetNrow() > 0);
        assert(GetNrow() == static_cast<int>(Ne.size()));
        for (int i = 0; i < this->GetNrow(); i++) {
            SetRowNe(i, Ne[i]);
        }
    }

    //! signal corruption of the parameters (matrices should recompute themselves)
    void CorruptCodonMatrices() {
        for (int j = 0; j < GetNcol(); j++) { CorruptColCodonMatrices(j); }
    }

    //! signal corruption for column (site) j
    void CorruptColCodonMatrices(int j) {
        assert(GetNrow() > 0);
        for (int i = 0; i < GetNrow(); i++) { CorruptMatrix(i, j); }
    }

    //! signal corruption for row (branch) i
    void CorruptRowCodonMatrices(int i) {
        assert(GetNcol() > 0);
        for (int j = 0; j < GetNcol(); j++) { CorruptMatrix(i, j); }
    }

    //! signal corruption for row (branch) i and column (site) j
    void CorruptMatrix(int i, int j) {
        (*this)(i, j).CorruptMatrix();
    }

    //! update only those matrices for which occupancy[i] != 0
    void UpdateCodonMatrices(const Selector<int> &occupancy) {
        for (int col = 0; col < GetNcol(); col++) {
            if (!occupancy.GetVal(col)) { (*this).CorruptColCodonMatrices(col); }
        }
    }

  private:
    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const Selector<std::vector<double>> *aafitnessarray;
    const std::vector<double> *popsizearray;
    std::vector<std::vector<AAMutSelOmegaCodonSubMatrix *>> matrixbidimarray;
};

/**
 * \brief An array of mutation-selection codon matrices (with Ne)
 *
 * The array takes a single nucleotide matrix, an array of amino-acid fitness
 * profiles and a population size (Ne) and construct
 * mutation-selection matrices accordingly.
 */

class AAMutSelNeCodonSubMatrixArray : public Array<SubMatrix>,
                                      public Array<AAMutSelOmegaCodonSubMatrix> {
  public:
    //! constructor with a nucleotide matrix, an array of amino-acid fitness
    //! profiles and a single Ne value (for all matrices)
    AAMutSelNeCodonSubMatrixArray(const CodonStateSpace *incodonstatespace,
        const SubMatrix *innucmatrix, const Selector<std::vector<double>> *inaafitnessarray,
        double inne)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          aafitnessarray(inaafitnessarray),
          ne(inne),
          matrixarray(inaafitnessarray->GetSize()) {
        Create();
    }

    ~AAMutSelNeCodonSubMatrixArray() { Delete(); }

    //! \brief set omega to new value
    //!
    //! should be called only when all matrices share same Ne parameter.
    //! makes an error (with exit) if this is not the case.
    void SetNe(double inne) { ne = inne; }

    //! return array size
    int GetSize() const { return aafitnessarray->GetSize(); }
    //! const access to matrix i
    const AAMutSelOmegaCodonSubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    AAMutSelOmegaCodonSubMatrix &operator[](int i) { return *matrixarray[i]; }

    //! const acess to nucleotide matrix
    const SubMatrix &GetNucMatrix() const { return *nucmatrix; }

    //! update all matrices
    void UpdateCodonMatrices() {
        for (int i = 0; i < GetSize(); i++) { UpdateCodonMatrices(i); }
    }

    //! update matrices of component i
    void UpdateCodonMatrices(int i) {
        (*this)[i].SetNe(ne);
        (*this)[i].CorruptMatrix();
    }

  private:
    void Create() {
        for (int i = 0; i < GetSize(); i++) {
            matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(
                codonstatespace, nucmatrix, aafitnessarray->GetVal(i), 1.0, ne);
        }
    }

    void Delete() {
        for (int i = 0; i < GetSize(); i++) { delete matrixarray[i]; }
    }

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const Selector<std::vector<double>> *aafitnessarray;
    double ne;
    std::vector<AAMutSelOmegaCodonSubMatrix *> matrixarray;
};