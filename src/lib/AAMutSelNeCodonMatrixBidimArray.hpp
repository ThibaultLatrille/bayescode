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
    //! return the number of columns (number of sites)
    int GetNcol() const override { return aafitnessarray->GetSize(); }

    //! const access to the matrix for row (branch) i and column (site) j
    const AAMutSelOmegaCodonSubMatrix &GetVal(int i, int j) const override;

    //! non const access to the matrix for row (branch) i and column (site) j
    AAMutSelOmegaCodonSubMatrix &operator()(int i, int j) override;

    //! Set ne for row (branch) i
    void SetRowNe(int i, double Ne);

    void SetNe(std::vector<double> const &Ne);

    //! signal corruption of the parameters (matrices should recompute themselves)
    void CorruptCodonMatrices();

    //! signal corruption for column (site) j
    void CorruptColCodonMatrices(int j);

    //! signal corruption for row (branch) i
    void CorruptRowCodonMatrices(int i);

    //! signal corruption for row (branch) i and column (site) j
    void CorruptMatrix(int i, int j);

    //! update only those matrices for which occupancy[i] != 0
    void UpdateCodonMatrices(const Selector<int> &occupancy);

  private:
    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const Selector<std::vector<double>> *aafitnessarray;
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
    void UpdateCodonMatrices();

    //! update matrices of component i
    void UpdateCodonMatrices(int i);

  private:
    void Create();

    void Delete();

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const Selector<std::vector<double>> *aafitnessarray;
    double ne;
    std::vector<AAMutSelOmegaCodonSubMatrix *> matrixarray;
};