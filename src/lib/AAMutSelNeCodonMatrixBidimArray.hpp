#pragma once

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "BidimArray.hpp"

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
        const std::vector<double> &pop_size_array);

    ~MutSelNeCodonMatrixBidimArray() override;

    //! return the number of rows (number of branches)
    int GetNrow() const override { return matrixbidimarray.size(); }
    //! return the number of columns (number of cat)
    int GetNcol() const override { return matrixbidimarray.begin()->size(); }

    //! const access to the matrix for row (branch) i and column (cat) j
    const AAMutSelOmegaCodonSubMatrix &GetVal(int i, int j) const override;

    //! non const access to the matrix for row (branch) i and column (cat) j
    AAMutSelOmegaCodonSubMatrix &operator()(int i, int j) override;

    //! Set ne for row (branch) i
    void UpdateRowNe(int i, double Ne);

    //! signal corruption when fitness and Ne have not changed
    void UpdateCodonMatricesNoFitnessRecomput();

    //! signal corruption for column (cat) j
    void UpdateColCodonMatrices(int j);

    //! update only those matrices for which occupancy[i] != 0
    void UpdateColCodonMatrices(const Selector<int> &occupancy);

    //! signal corruption for row (branch) i and column (cat) j
    void UpdateMatrix(int i, int j);

  private:
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
        double inne);

    ~AAMutSelNeCodonSubMatrixArray() override { Delete(); }

    //! return array size
    int GetSize() const override { return matrixarray.size(); }
    //! const access to matrix i
    const AAMutSelOmegaCodonSubMatrix &GetVal(int i) const override { return *matrixarray[i]; }
    //! non-const access to matrix i
    AAMutSelOmegaCodonSubMatrix &operator[](int i) override { return *matrixarray[i]; }

    //! \brief set Ne to new value
    //!
    //! should be called only when all matrices share same Ne parameter.
    //! makes an error (with exit) if this is not the case.
    void UpdateNe(double ne) {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].UpdateNe(ne); }
    }

    //! update all matrices
    void UpdateCodonMatricesNoFitnessRecomput() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].CorruptMatrixNoFitnessRecomput(); }
    };

    //! update matrices of component i
    void UpdateCodonMatrices(int i) { (*this)[i].CorruptMatrix(); };

    void UpdateCodonMatrices(const Selector<int> &occupancy) {
        for (int i = 0; i < GetSize(); i++) {
            if (!occupancy.GetVal(i)) { this->UpdateCodonMatrices(i); }
        }
    }

  private:
    void Delete();

    std::vector<AAMutSelOmegaCodonSubMatrix *> matrixarray;
};