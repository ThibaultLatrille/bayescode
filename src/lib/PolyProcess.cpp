#include "PolyProcess.hpp"
#include <numeric>

using namespace std;

PolyProcess::PolyProcess(CodonStateSpace &instatespace, PolyData &indata,
    PoissonRandomField &inpoissonrandomfield, Selector<vector<double>> const &insiteaafitnessarray,
    GTRSubMatrix const &innucmatrix, ScaledMutationRate const &intheta)
    : polydata{indata},
      poissonrandomfield{inpoissonrandomfield},
      statespace{instatespace},
      nucmatrix{innucmatrix},
      siteaafitnessarray{insiteaafitnessarray},
      theta{intheta} {}

double PolyProcess::GetProb(int taxon, int site, int anc_state) {
    unsigned sample_size = polydata.GetSampleSize(taxon, site);
    unsigned anc_occurence = polydata.GetCount(taxon, site, anc_state);

    if (anc_occurence == sample_size) {
        // If the ancestral allele is monomorphic
        return poissonrandomfield.GetProb(anc_state, anc_state, sample_size, sample_size,
            siteaafitnessarray.GetVal(site), nucmatrix, theta.GetTheta(taxon));
    } else {
        // If the ancestral allele is not monomorphic

        for (auto der_state : statespace.GetNeighbors(anc_state)) {
            unsigned der_occurence = polydata.GetCount(taxon, site, der_state);
            if (der_occurence + anc_occurence == sample_size) {
                return poissonrandomfield.GetProb(anc_state, der_state, der_occurence, sample_size,
                    siteaafitnessarray.GetVal(site), nucmatrix, theta.GetTheta(taxon));
            }
        }
    }
    return 0.0;
}

double PolyProcess::GetLogProb(int taxon, int site, int anc_state) {
    double proba = GetProb(taxon, site, anc_state);
    if (proba > 0) {
        return log(proba);
    } else {
        return -numeric_limits<double>::infinity();
    }
}

tuple<int, int, unsigned, unsigned> PolyProcess::GetDerivedTuple(
    int taxon, int site, int anc_state) {
    unsigned sample_size = polydata.GetSampleSize(taxon, site);
    unsigned anc_occurence = polydata.GetCount(taxon, site, anc_state);

    if (anc_occurence == sample_size) {
        // If the ancestral allele is monomorphic
        return make_tuple(anc_state, anc_state, sample_size, sample_size);
    } else {
        // If the ancestral allele is not monomorphic

        for (auto der_state : statespace.GetNeighbors(anc_state)) {
            unsigned der_occurence = polydata.GetCount(taxon, site, der_state);
            if (der_occurence + anc_occurence == sample_size) {
                return make_tuple(anc_state, der_state, der_occurence, sample_size);
            }
        }
    }
    cerr << "GetDerivedTuple should have returned" << endl;
    return make_tuple(-1, -1, 0, 0);
}