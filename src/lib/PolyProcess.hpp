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
    PolyProcess(CodonStateSpace &instatespace, PolyData &indata,
        PoissonRandomField &inpoissonrandomfield,
        Selector<std::vector<double>> const &insiteaafitnessarray, GTRSubMatrix const &innucmatrix,
        ScaledMutationRate const &intheta);

    ~PolyProcess() /*override*/ = default;

    //! \brief give the likelihood of the data.
    //! This method can't be declared const because evaluation of the
    //! likelihood might need to update the pre-computed values
    double GetProb(int taxon, int site, int anc_state);

    //! Log likelihood
    double GetLogProb(int taxon, int site, int anc_state);

    std::tuple<int, int, unsigned, unsigned> GetDerivedTuple(int taxon, int site, int anc_state);

    CodonStateSpace *GetCodonStateSpace() { return &statespace; };
    PolyData *GetPolyData() { return &polydata; };

  private:
    PolyData &polydata;
    PoissonRandomField &poissonrandomfield;
    CodonStateSpace &statespace;

    GTRSubMatrix const &nucmatrix;
    Selector<std::vector<double>> const &siteaafitnessarray;
    ScaledMutationRate const &theta;
};
