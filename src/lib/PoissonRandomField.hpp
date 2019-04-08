#pragma once

#include <deque>
#include <set>
#include <vector>
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

/**
 * \brief PolyProcess is used by PhyloProcess to compute the likelihood the data, in the case
 * where polymorphism dataset are available.
 */

class PoissonRandomField {
  public:
    //! \brief Constructor: takes a (const pointer to a) PolyData
    PoissonRandomField(std::set<unsigned> const &sample_size_set, CodonStateSpace const &instatespace,
        unsigned precision);

    ~PoissonRandomField() /*override*/ = default;

    //! \brief give the likelihood of the data.
    //! This method can't be declared const because evaluation of the
    //! likelihood might need to update the pre-computed values
    double GetProb(int anc_state, int der_state, unsigned der_occurence, unsigned sample_size,
        const std::vector<double> &aafitnessarray, const GTRSubMatrix &nucmatrix,
        const double &theta) const;


  private:
    CodonStateSpace const &statespace;

    unsigned precision{10};
    double grid_s_step{0};

    std::map<unsigned, std::deque<std::pair<double, std::vector<double>>>> mutable ComputedProb;
    std::map<unsigned, std::vector<unsigned long long>> mutable ComputedBinom;

    //! Use the pre-computed values, and update the data if interpolation is not possible
    double InterpolateProba(int anc_state, int der_state, unsigned der_occurence,
        unsigned sample_size, const std::vector<double> &aafitnessarray,
        const GTRSubMatrix &nucmatrix) const;

    //! Update the precomputed values
    void UpdateComputed(unsigned sample_size, double s) const;

    //! Let i be the number of copies of the derived allele, in a sample of size n.
    //! Return the vector (for all i from 1 to n) of expected time for which we
    //! observe i copies in the sample of size n.
    //! The derived allele has a selection coefficient (s).
    std::vector<double> ExpectedTimeObsVector(unsigned n, double s) const;
};

//! Return a string from a vector of T (double, int, ..),
//! values are also separated by a char (by default a whitespace)
template <class T>
std::string join(std::vector<T> const &v, char sep = ' ');

//! Return the vector of binomial coefficient from 0 to n (included)
std::vector<unsigned long long> BinomialCoefficientArray(unsigned n);

//! Return a exponent b (a^b) where a and b are unsigned values
unsigned PowUnsigned(unsigned a, unsigned b);

//! Return the sum of a vector of double
double SumVector(std::vector<double> const &v);
