#pragma once

#include "BranchArray.hpp"
#include "CodonSubMatrix.hpp"
#include "SubMatrix.hpp"

/**
 * \brief A BranchArray of MGOmegaCodonSubMatrix
 *
 * The constructor takes as arguments a single nucleotide matrix and a
 * BranchSelector<double> which returns the value of omega for each site. It
 * then constructs an array of MGOmegaCodonSubMatrix of same size as the omega
 * array, and such that the matrix on branch j takes the jth value of the omega
 * array as its dN/dS.
 *
 * This array derives from BranchArray and not SimpleBranchArray, because it is
 * implemented as a vector<MGOmegaCodonSubMatrix*>.
 */

class MGOmegaCodonSubMatrixBranchArray : public BranchArray<SubMatrix>,
                                         public BranchArray<MGOmegaCodonSubMatrix> {
  public:
    //! constructor parameterized by a codon state space, a single nucleotide
    //! matrix and an array (in fact, a BranchSelector) of omega's
    MGOmegaCodonSubMatrixBranchArray(const CodonStateSpace *incodonstatespace,
        const SubMatrix *innucmatrix, const BranchSelector<double> *inomegaarray)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          omegaarray(inomegaarray),
          matrixarray(inomegaarray->GetNbranch()) {
        Create();
    }

    ~MGOmegaCodonSubMatrixBranchArray() { Delete(); }

    //! return array size
    const Tree &GetTree() const { return omegaarray->GetTree(); }
    //! const access to matrix i
    int GetNbranch() const { return omegaarray->GetNbranch(); }
    const MGOmegaCodonSubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    MGOmegaCodonSubMatrix &operator[](int i) { return *matrixarray[i]; }

    //! const access to underlying nucleotide matrix
    const SubMatrix &GetNucMatrix() const { return *nucmatrix; }

    //! update all matrices
    void UpdateCodonMatrices() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].UpdateOmega(omegaarray->GetVal(i));
            (*this)[i].CorruptMatrix();
        }
    }

  private:
    void Create() {
        for (int i = 0; i < GetNbranch(); i++) {
            matrixarray[i] =
                new MGOmegaCodonSubMatrix(codonstatespace, nucmatrix, omegaarray->GetVal(i));
        }
    }

    void Delete() {
        for (int i = 0; i < GetNbranch(); i++) { delete matrixarray[i]; }
    }

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const BranchSelector<double> *omegaarray;
    std::vector<MGOmegaCodonSubMatrix *> matrixarray;
};