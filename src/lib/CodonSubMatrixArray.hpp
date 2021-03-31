#pragma once

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "Array.hpp"
#include "CodonSubMatrix.hpp"
#include "SubMatrix.hpp"

/**
 * \brief An Array of MGOmegaCodonSubMatrix
 *
 * The constructor takes as arguments a single nucleotide matrix and a
 * Selector<double> which returns the value of omega for each site. It then
 * constructs an array of MGOmegaCodonSubMatrix of same size as the omega array,
 * and such that the ith matrix takes the ith value of the omega array as its
 * dN/dS.
 *
 * This array derives from Array and not SimpleArray, because it is implemented
 * as a vector<MGOmegaCodonSubMatrix*>.
 */

class MGOmegaCodonSubMatrixArray : public Array<SubMatrix>, public Array<MGOmegaCodonSubMatrix> {
  public:
    //! constructor parameterized by a codon state space, a single nucleotide
    //! matrix and an array (in fact, a Selector) of omega's
    MGOmegaCodonSubMatrixArray(const CodonStateSpace *incodonstatespace,
        const SubMatrix *innucmatrix, const Selector<double> *inomegaarray)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          omegaarray(inomegaarray),
          matrixarray(inomegaarray->GetSize()) {
        Create();
    }

    ~MGOmegaCodonSubMatrixArray() { Delete(); }

    //! return array size
    int GetSize() const { return omegaarray->GetSize(); }
    //! const access to matrix i
    const MGOmegaCodonSubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    MGOmegaCodonSubMatrix &operator[](int i) { return *matrixarray[i]; }

    //! const access to underlying nucleotide matrix
    const SubMatrix &GetNucMatrix() const { return *nucmatrix; }

    //! update all matrices
    void UpdateCodonMatrices() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].UpdateOmega(omegaarray->GetVal(i));
            (*this)[i].CorruptMatrix();
        }
    }

    //! update only those matrices for which occupancy[i] != 0
    void UpdateCodonMatrices(const Selector<int> &occupancy) {
        if (occupancy.GetSize() != GetSize()) {
            std::cerr << "error in UpdateCodonMatrices: occupancy vector size does not "
                         "match array size\n";
            exit(1);
        }
        for (int i = 0; i < GetSize(); i++) {
            if (occupancy.GetVal(i)) {
                (*this)[i].UpdateOmega(omegaarray->GetVal(i));
                (*this)[i].CorruptMatrix();
            }
        }
    }

  private:
    void Create() {
        for (int i = 0; i < GetSize(); i++) {
            matrixarray[i] =
                new MGOmegaCodonSubMatrix(codonstatespace, nucmatrix, omegaarray->GetVal(i));
        }
    }

    void Delete() {
        for (int i = 0; i < GetSize(); i++) { delete matrixarray[i]; }
    }

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const Selector<double> *omegaarray;
    std::vector<MGOmegaCodonSubMatrix *> matrixarray;
};

/**
 * \brief An array of mutation-selection codon matrices (with omega)
 *
 * The array takes a single nucleotide matrix, an array of amino-acid fitness
 * profiles and either a single omega or an array of omega values and construct
 * mutation-selection matrices accordingly.
 */

class AAMutSelOmegaCodonSubMatrixArray : public Array<SubMatrix>,
                                         public Array<AAMutSelOmegaCodonSubMatrix> {
  public:
    //! constructor with a nucleotide matrix, an array of amino-acid fitness
    //! profiles and a single omega value (for all matrices)
    AAMutSelOmegaCodonSubMatrixArray(const CodonStateSpace *incodonstatespace,
        const SubMatrix *innucmatrix, const Selector<std::vector<double>> *inaafitnessarray,
        double inomega)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          aafitnessarray(inaafitnessarray),
          omega(inomega),
          omegaarray(0),
          matrixarray(inaafitnessarray->GetSize()) {
        Create();
    }

    //! constructor with a nucleotide matrix, an array of amino-acid fitness
    //! profiles and an array of omega value (one for each entry of the matrix
    //! array)
    AAMutSelOmegaCodonSubMatrixArray(const CodonStateSpace *incodonstatespace,
        const SubMatrix *innucmatrix, const Selector<std::vector<double>> *inaafitnessarray,
        const Selector<double> *inomegaarray)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          aafitnessarray(inaafitnessarray),
          omegaarray(inomegaarray),
          matrixarray(inomegaarray->GetSize()) {
        if (aafitnessarray->GetSize() != omegaarray->GetSize()) {
            std::cerr << "error in constructor of AAMutSelOmegaCodonSubMatrixArray: "
                         "arrays of aafitness and omega values should be of same size\n";
            exit(1);
        }
        Create();
    }

    ~AAMutSelOmegaCodonSubMatrixArray() { Delete(); }

    //! \brief set omega to new value
    //!
    //! should be called only when all matrices share same omega parameter.
    //! makes an error (with exit) if this is not the case.
    void SetOmega(double inomega) {
        if (omegaarray) {
            std::cerr << "error in AAMutSelOmegaCodonSubMatrixArray::UpdateOmega\n";
            exit(1);
        }
        omega = inomega;
    }

    //! return array size
    int GetSize() const { return aafitnessarray->GetSize(); }
    //! const access to matrix i
    const AAMutSelOmegaCodonSubMatrix &GetVal(int i) const { return *matrixarray[i]; }
    //! non-const access to matrix i
    AAMutSelOmegaCodonSubMatrix &operator[](int i) { return *matrixarray[i]; }

    //! const acess to nucleotide matrix
    const SubMatrix &GetNucMatrix() const { return *nucmatrix; }

    //! update all matrices
    void UpdateCodonMatrices(bool fitness_recomput = true) {
        for (int i = 0; i < GetSize(); i++) {
            if (omegaarray) {
                (*this)[i].UpdateOmega(omegaarray->GetVal(i));
            } else {
                (*this)[i].UpdateOmega(omega);
            }
            if (fitness_recomput) {
                (*this)[i].CorruptMatrix();
            } else {
                (*this)[i].CorruptMatrixNoFitnessRecomput();
            }
        }
    }

    //! update only those matrices for which occupancy[i] != 0
    void UpdateCodonMatrices(const Selector<int> &occupancy) {
        if (omegaarray) {
            for (int i = 0; i < GetSize(); i++) {
                if (!occupancy.GetVal(i)) {
                    (*this)[i].UpdateOmega(omegaarray->GetVal(i));
                    (*this)[i].CorruptMatrix();
                }
            }
        } else {
            for (int i = 0; i < GetSize(); i++) {
                if (!occupancy.GetVal(i)) {
                    (*this)[i].UpdateOmega(omega);
                    (*this)[i].CorruptMatrix();
                }
            }
        }
    }

  private:
    void Create() {
        for (int i = 0; i < GetSize(); i++) {
            if (omegaarray) {
                matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(codonstatespace, nucmatrix,
                    aafitnessarray->GetVal(i), omegaarray->GetVal(i), 1.0);
            } else {
                matrixarray[i] = new AAMutSelOmegaCodonSubMatrix(
                    codonstatespace, nucmatrix, aafitnessarray->GetVal(i), omega, 1.0);
            }
        }
    }

    void Delete() {
        for (int i = 0; i < GetSize(); i++) { delete matrixarray[i]; }
    }

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const Selector<std::vector<double>> *aafitnessarray;
    double omega;
    const Selector<double> *omegaarray;
    std::vector<AAMutSelOmegaCodonSubMatrix *> matrixarray;
};