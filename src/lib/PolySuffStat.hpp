#pragma once

#include <map>
#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchArray.hpp"
#include "NodeArray.hpp"
#include "PhyloProcess.hpp"
#include "PoissonRandomField.hpp"
#include "SubMatrix.hpp"
#include "SuffStat.hpp"

/**
 * \brief A Poisson-Random-Field sufficient statistic
 */

class PolySuffStat : public SuffStat {
  public:
    PolySuffStat() {}

    ~PolySuffStat() {}

    //! set suff stats to 0
    void Clear() { polycount.clear(); }

    void IncrementPolyCount(std::tuple<int, int, unsigned, unsigned> poly_tuple) {
        polycount[poly_tuple]++;
    }

    void AddPairCount(std::tuple<int, int, unsigned, unsigned> poly_tuple, int in) {
        polycount[poly_tuple] += in;
    }

    //! add Poly sufficient statistics from PhyloProcess (site-homogeneous case)
    void AddSuffStat(const PhyloProcess &process) { process.AddPolySuffStat(*this); }

    void Add(const PolySuffStat &suffstat) {
        for (auto const &i : suffstat.GetPairCountMap()) { AddPairCount(i.first, i.second); }
    }

    PolySuffStat &operator+=(const PolySuffStat &from) {
        Add(from);
        return *this;
    }

    int GetPairCount(std::tuple<int, int, unsigned, unsigned> poly_tuple) const {
        auto const i = polycount.find(poly_tuple);
        if (i == polycount.end()) { return 0; }
        return i->second;
    }

    //! return log p(S | Q) as a function of the fitness vector, the nucmatrix and theta
    double GetLogProb(PoissonRandomField &poissonrandomfield,
        const std::vector<double> &aafitnessarray, const GTRSubMatrix &nucmatrix,
        const double &theta) const {
        double total = 0;
        for (auto const &i : polycount) {
            double proba = poissonrandomfield.GetProb(std::get<0>(i.first), std::get<1>(i.first),
                std::get<2>(i.first), std::get<3>(i.first), aafitnessarray, nucmatrix, theta);
            if (proba > 0) {
                total += i.second * log(proba);
            } else {
                total = -std::numeric_limits<double>::infinity();
            }
        }
        return total;
    }

    //! const access to the ordered map giving the pair count stat (sparse data
    //! structure)
    const std::map<std::tuple<int, int, unsigned, unsigned>, int> &GetPairCountMap() const {
        return polycount;
    }

  private:
    std::map<std::tuple<int, int, unsigned, unsigned>, int> polycount;
};

/**
 * \brief An array of Poisson-Random-Field sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each item (each site, or each component of a mixture) has a different rate
 * matrix Q_i
 */
class PolySuffStatArray : public SimpleArray<PolySuffStat> {
  public:
    PolySuffStatArray(int insize) : SimpleArray<PolySuffStat>(insize) {}

    ~PolySuffStatArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Clear(); }
    }

    //! add Poly sufficient statistics from PhyloProcess (site-heterogeneous case)
    void AddSuffStat(const PhyloProcess &process) { process.AddPolySuffStat(*this); }

    //! return total log prob (summed over all items), given an array of fitness vector, the
    //! nucmatrix and theta
    double GetLogProb(PoissonRandomField &poissonrandomfield,
        const Selector<std::vector<double>> &siteaafitnessarray, const GTRSubMatrix &nucmatrix,
        const double &theta) const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetVal(i).GetLogProb(
                poissonrandomfield, siteaafitnessarray.GetVal(i), nucmatrix, theta);
        }
        return total;
    }

    //! \brief add suffstatarray given as argument to this array based on the
    //! allocations provided as the second argument (mixture models)
    //!
    void Add(const Selector<PolySuffStat> &suffstatarray, const Selector<int> &alloc) {
        for (int i = 0; i < suffstatarray.GetSize(); i++) {
            (*this)[alloc.GetVal(i)] += suffstatarray.GetVal(i);
        }
    }
};

/**
 * \brief A 2-dimensional-array of Poisson-Random-Field sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each item (each site, or each component of a mixture and each branch or each component of the
 * branch allocation) has a different rate matrix Q_(i,j)
 */
class PolySuffStatBidimArray : public SimpleBidimArray<PolySuffStat> {
  public:
    PolySuffStatBidimArray(int inrow, int incol)
        : SimpleBidimArray<PolySuffStat>(inrow, incol, PolySuffStat()) {}

    ~PolySuffStatBidimArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int row = 0; row < GetNrow(); row++) {
            for (int col = 0; col < GetNcol(); col++) { (*this)(row, col).Clear(); }
        }
    }

    //! add Poly sufficient statistics from PhyloProcess (site-heterogeneous and
    //! branch-heterogeneous case)
    void AddSuffStat(const PhyloProcess &process) { process.AddPolySuffStat(*this); }

    //! return total log prob (summed over all items), given an array of fitness vector, the
    //! nucmatrix and theta per taxon
    double GetLogProb(PoissonRandomField &poissonrandomfield,
        const Selector<std::vector<double>> &siteaafitnessarray, const GTRSubMatrix &nucmatrix,
        const ScaledMutationRate &theta) const {
        double total = 0;
        for (int row = 0; row < GetNrow(); row++) {
            double d_theta = theta.GetTheta(row);
            for (int col = 0; col < GetNcol(); col++) {
                total += GetVal(row, col).GetLogProb(poissonrandomfield,
                    siteaafitnessarray.GetVal(col), nucmatrix, d_theta);
            };
        }
        return total;
    }

    //! \brief add suffstatarray given as argument to this array based on the
    //! allocations provided as the second argument (mixture models)
    //!
    void Add(
        const BidimSelector<PolySuffStat> &suffstatbidimarray, const Selector<int> &col_alloc) {
        for (int row = 0; row < suffstatbidimarray.GetNrow(); row++) {
            for (int col = 0; col < suffstatbidimarray.GetNcol(); col++) {
                (*this)(row, col_alloc.GetVal(col)) += suffstatbidimarray.GetVal(row, col);
            }
        }
    }
};