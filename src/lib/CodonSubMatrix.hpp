#pragma once

#include "CodonStateSpace.hpp"
#include "SubMatrix.hpp"

const double omegamin = 1e-5;

/**
 * \brief A general class representing all codon matrices
 *
 * Provides a few generic methods for dealing with codon-specific issues (syn,
 * non syn, etc). However, still needs to be specialized into actual classes
 * (the two most important being MGOmegaCodonSubMatrix and
 * AAMutSelOmegaCodonSubMatrix)
 */

class CodonSubMatrix : public virtual SubMatrix {
  public:
    //! constructor parameterized by codon state space (itself specifying the
    //! genetic code)
    CodonSubMatrix(const CodonStateSpace *instatespace, bool innormalise)
        : SubMatrix(instatespace->GetNstate(), innormalise), statespace(instatespace) {}

    const CodonStateSpace *GetCodonStateSpace() const { return statespace; }

    //! see CodonStateSpace::Synonymous
    bool Synonymous(int codon1, int codon2) const { return statespace->Synonymous(codon1, codon2); }
    //! see CodonStateSpace::GetCodonPosition
    int GetCodonPosition(int pos, int codon) const {
        return statespace->GetCodonPosition(pos, codon);
    }
    //! see CodonStateSpace::GetDifferingPosition
    int GetDifferingPosition(int codon1, int codon2) const {
        return statespace->GetDifferingPosition(codon1, codon2);
    }

  protected:
    const CodonStateSpace *statespace;
};

/**
 * \brief A generic class for codon matrices based on an underlying nucleotide
 * rate matrix
 *
 * Most codon matrices rely on a mutation process at the level of nucleotides.
 * NucCodonSubMatrix takes a const pointer to a nucleotide substitution matrix
 * representing the nucleotide 4x4 mutation process.
 * Still an abstract class.
 */
class NucCodonSubMatrix : public virtual CodonSubMatrix {
  public:
    NucCodonSubMatrix(
        const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix, bool innormalise)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise) {
        SetNucMatrix(inNucMatrix);
    }

    const SubMatrix *GetNucMatrix() const { return NucMatrix; }

  protected:
    void SetNucMatrix(const SubMatrix *inmatrix) {
        NucMatrix = inmatrix;
        if (NucMatrix->GetNstate() != Nnuc) {
            std::cerr << "error in CodonSubMatrix: underyling mutation process "
                         "should be a 4x4 matrix\n";
            throw;
        }
    }

    const SubMatrix *NucMatrix;
};

/**
 * \brief A generic class for codon matrices with an omega parameter, acting as
 * a multiplier in front of all rates between non-synonynmous codons
 */

class OmegaCodonSubMatrix : public virtual CodonSubMatrix {
  public:
    OmegaCodonSubMatrix(const CodonStateSpace *instatespace, double inomega, bool innormalise)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, normalise),
          omega(inomega) {}

    double GetOmega() const { return omega + omegamin; }
    void UpdateOmega(double inomega) {
        omega = inomega;
        CorruptMatrix();
    }

  protected:
    double omega;
};

/**
 * \brief A Muse and Gaut codon substitution process (Muse and Gaut, 1994).
 *
 * The simplest codon model based on a pure nucleotide mutation process (with
 * stops excluded). Without omega.
 */

class MGCodonSubMatrix : public NucCodonSubMatrix {
  public:
    MGCodonSubMatrix(
        const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise) {}

    void CorruptMatrix() /*override*/ { SubMatrix::CorruptMatrix(); }

  protected:
    void ComputeArray(int i) const /*override*/;
    void ComputeStationary() const /*override*/;
};

/**
 * \brief A Muse and Gaut codon substitution process with an omgea = dN/dS
 * parameter
 */

class MGOmegaCodonSubMatrix : public MGCodonSubMatrix, public OmegaCodonSubMatrix {
  public:
    MGOmegaCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
        double inomega, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          OmegaCodonSubMatrix(instatespace, inomega, innormalise) {}

  protected:
    void ComputeArray(int i) const /*override*/;
};