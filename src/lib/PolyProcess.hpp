#pragma once

#include <deque>
#include "Array.hpp"
#include "GTRSubMatrix.hpp"
#include "PoissonRandomField.hpp"
#include "PolyData.hpp"
#include "ScaledMutationRate.hpp"

/**
 * \brief PolyProcess is used by PhyloProcess to compute the likelihood the data, in the case
 * where polymorphism dataset are available.
 */

class PolyProcess {
  public:
    //! \brief Constructor: takes a (const pointer to a) PolyData
    PolyProcess(CodonStateSpace const &instatespace, PolyData const &indata,
        PoissonRandomField const &inpoissonrandomfield,
        MixtureSelector<std::vector<double>> const &insiteaafitnessarray,
        GTRSubMatrix const &innucmatrix, ScaledMutationRate const &intheta);

    ~PolyProcess() /*override*/ = default;

    //! \brief give the likelihood of the data.
    //! This method can't be declared const because evaluation of the
    //! likelihood might need to update the pre-computed values
    double GetProb(int taxon, int site, int anc_state) const;

    //! Log likelihood
    double GetLogProb(int taxon, int site, int anc_state) const;

    std::tuple<int, int, unsigned, unsigned> GetDerivedTuple(int taxon, int site, int anc_state) const;

  private:
    PolyData const &polydata;
    PoissonRandomField const &poissonrandomfield;
    CodonStateSpace const &statespace;

    GTRSubMatrix const &nucmatrix;
    MixtureSelector<std::vector<double>> const &siteaafitnessarray;
    ScaledMutationRate const &theta;
};
