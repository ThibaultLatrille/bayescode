#pragma once

#include "CodonSubMatrix.hpp"

/**
 * \brief A mutation-selection codon substitution process.
 *
 * This codon substitution process describes the evolution of a coding position,
 * under a constant fitness landscape over the 20 amino-acids. The model is
 * parameterized by a nucleotide matrix, specifying the mutation process, and a
 * vector of 20 scaled fitness parameters (summing to 1), for the 20
 * amino-acids. The process also takes a real parameter, omega, which acts as a
 * multiplier in front of all non-synonymous substitutions. The standard
 * mutation-selection process is obtained by setting omega=1. Letting omega be
 * different from 1 was explored in Rodrigue and Lartillot, 2017, MBE (detecting
 * deviations from expected non-syn rate under the standard mut-sel model).
 */

class AAMutSelOmegaCodonSubMatrix : public virtual NucCodonSubMatrix,
                                    public virtual OmegaCodonSubMatrix {
  public:
    //! constructor, parameterized by a codon state space (genetic code), a
    //! nucleotide mutation matrix, a 20-vector of amino-acid fitnesss, and a
    //! positive real parameter omega (=1 in the standard model).
    AAMutSelOmegaCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
        const std::vector<double> &inaa, double inomega, double inNe, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          OmegaCodonSubMatrix(instatespace, inomega, innormalise),
          aa(inaa),
          Ne(inNe) {}

    //! const access (by reference) to amino-acid fitness vector
    const std::vector<double> &GetAAFitnessProfile() const { return aa; }

    //! \brief access by copy to fitness of a given amino-acid
    //!
    //! Note: to avoid numerical errors, this function returns aa[a] + 1e-8.
    double GetFitness(int a) const { return exp(Ne * log(aa[a])) + 1e-8; }
    double GetLogFitness(int a) const { return log(GetFitness(a)); }

    std::tuple<double, double> GetFlowDNDS() const;
    double GetPredictedDNDS() const;

    void SetNe(double inNe) {
        Ne = inNe;
        CorruptMatrix();
    }

  protected:

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;


    // data members

    const std::vector<double> &aa;
    double Ne;
};