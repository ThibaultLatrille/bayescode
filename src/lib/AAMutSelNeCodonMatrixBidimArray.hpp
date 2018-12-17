#pragma once

#include <cassert>
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "Array.hpp"
#include "BidimArray.hpp"
#include "IIDDirichlet.hpp"

/**
 * \brief A BidimArray of AAMutSelOmegaCodonSubMatrix
 *
 * The constructor takes as arguments a single nucleotide matrix and a bidim array of amino-acid
 * fitness profiles (across sites and branches). It then constructs an array of
 * AAMutSelOmegaCodonSubMatrix of same size as the array of fitness profiles.
 */

class DatedMutSelCodonMatrixBidimArray : public BidimArray<SubMatrix>,
                                              public BidimArray<AAMutSelOmegaCodonSubMatrix> {
  public:
    //! constructor parameterized by a bidim array of fitness profiles, a codon state space and a
    //! single nucleotide matrix.
    DatedMutSelCodonMatrixBidimArray(const CodonStateSpace* incodonstatespace,
        const SubMatrix* innucmatrix, const Selector<std::vector<double>>* infitnessarray,
        const std::vector<double>* inNe, double inomega)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          aafitnessarray(infitnessarray),
          condition_ne(inNe),
          matrixbidimarray(inNe->size(),
              std::vector<AAMutSelOmegaCodonSubMatrix*>(infitnessarray->GetSize(), nullptr)) {
        std::cout << GetNrow() << "\t" << GetNcol() << "\n";
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                matrixbidimarray[i][j] = new AAMutSelOmegaCodonSubMatrix(codonstatespace, nucmatrix,
                    infitnessarray->GetVal(j), 1.0, condition_ne->at(i));
            }
        }
        assert(static_cast<int>(matrixbidimarray.size()) == GetNrow());
    }

    ~DatedMutSelCodonMatrixBidimArray() override {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) { delete matrixbidimarray[i][j]; }
        }
    }

    //! return the number of rows (number of branches)
    int GetNrow() const override { return condition_ne->size(); }
    //! return the number of columns (number of sites)
    int GetNcol() const override { return aafitnessarray->GetSize(); }

    //! const access to the matrix for row (branch) i and column (site) j
    const AAMutSelOmegaCodonSubMatrix& GetVal(int i, int j) const override {
        assert(0 <= i and i < GetNrow());
        assert(0 <= j and j < GetNcol());
        return *matrixbidimarray.at(i).at(j);
    }

    //! non const access to the matrix for row (branch) i and column (site) j
    AAMutSelOmegaCodonSubMatrix& operator()(int i, int j) override {
        assert(0 <= i and i < GetNrow());
        assert(0 <= j and j < GetNcol());
        return *matrixbidimarray.at(i).at(j);
    }

    void SetNe(std::vector<double>& Ne) {
        assert(GetNcol() > 0);
        assert(GetNrow() > 0);
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) { (*this)(i, j).SetNe(Ne[i]); }
        }
    }

    //! signal corruption of the parameters (matrices should recompute themselves)
    void UpdateCodonMatrices() {
        for (int j = 0; j < GetNcol(); j++) { UpdateCodonMatrices(j); }
    }

    //! signal corruption for column (site) j
    void UpdateCodonMatrices(int j) {
        assert(0 <= j and j < GetNcol());
        assert(GetNcol() > 0);
        assert(GetNrow() > 0);
        for (int i = 0; i < GetNrow(); i++) { CorruptMatrix(i, j); }
    }

    //! signal corruption for row (branch) i and column (site) j
    void CorruptMatrix(int i, int j) {
        assert(0 <= i and i < GetNrow());
        assert(GetNcol() > 0);
        assert(GetNrow() > 0);
        (*this)(i, j).CorruptMatrix();
    }

    //! signal corruption for column (site) j, and only for those rows (branches) that are flagged
    void UpdateCodonMatrices(int j, const std::vector<int>& flag) {
        assert(GetNcol() > 0);
        assert(GetNrow() > 0);
        assert(0 <= j and j < GetNcol());
        for (int i = 0; i < GetNrow(); i++) {
            if (flag[i]) { CorruptMatrix(i, j); }
        }
    }

    //! update only those matrices for which occupancy[i] != 0
    void UpdateCodonMatrices(const Selector<int>& occupancy) {
        for (int col = 0; col < GetNcol(); col++) {
            if (!occupancy.GetVal(col)) { (*this).UpdateCodonMatrices(col); }
        }
    }

  private:
    const CodonStateSpace* codonstatespace;
    const SubMatrix* nucmatrix;
    const Selector<std::vector<double>>* aafitnessarray;
    const std::vector<double>* condition_ne;
    std::vector<std::vector<AAMutSelOmegaCodonSubMatrix*>> matrixbidimarray;
};
