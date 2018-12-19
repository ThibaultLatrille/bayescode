#pragma once

#include "Array.hpp"
#include "BranchArray.hpp"
#include "PoissonSuffStat.hpp"
#include "Random.hpp"
#include "components/traits.hpp"

/**
 * \brief An array of IID gamma random variables
 */

class IIDGamma : public SimpleArray<double> {
public:
    //! constructor specifies array size and initial shape and scale parameters
    IIDGamma(int insize, const double &inmean, const double &ininvshape)
        : SimpleArray<double>(insize), mean(inmean), invshape(ininvshape) {
        Sample();
    }

    ~IIDGamma() {}

    //! return shape parameter
    double GetShape() const { return 1.0 / invshape; }

    //! return scale parameter
    double GetScale() const { return (1.0 / invshape) / mean; }

    //! sample all entries, given current shape and scale params
    void Sample() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] = Random::GammaSample(GetShape(), GetScale());
        }
    }

    //! resample all entries, given current shape and scale parameters and given
    //! an array of Poisson sufficient statistics of same size
    void GibbsResample(const Selector<PoissonSuffStat> &suffstatarray) {
        for (int i = 0; i < GetSize(); i++) {
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);
            (*this)[i] = Random::GammaSample(
                GetShape() + suffstat.GetCount(), GetScale() + suffstat.GetBeta());
        }
    }

    //! return total log prob (summed over the array) given current shape and
    //! scale params
    double GetLogProb() const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! return log prob of one specific entry
    double GetLogProb(int index) const {
        return Random::logGammaDensity(GetVal(index), GetShape(), GetScale());
    }

    //! \brief given a Poisson suffstat S, calculates postprob[i] propto weight[i]
    //! * p(S | (*this)[i]), for i=0..GetSize()-1.
    //!
    //! interpretation: postprob[i] = posterior probability that the data
    //! summarized by S have been produced by a process of rate (*this)[i], for
    //! i=0..GetSize()-1
    void GetAllocPostProb(const PoissonSuffStat &suffstat, const std::vector<double> &weight,
                          std::vector<double> &postprob) const {
        double max = 0;
        for (int i = 0; i < GetSize(); i++) {
            double tmp = suffstat.GetLogProb(GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) { max = tmp; }
        }

        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            postprob[i] = weight[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < GetSize(); i++) { postprob[i] /= total; }
    }

    //! resample all entries for which occupancy[i] == 0 from the prior (from a
    //! Gamma(shape,scale))
    void PriorResample(const Selector<int> &occupancy) {
        for (int i = 0; i < GetSize(); i++) {
            if (!occupancy.GetVal(i)) { (*this)[i] = Random::GammaSample(GetShape(), GetScale()); }
        }
    }

    //! resample all entries for which poswarray[i] == 0 from the prior (from a
    //! Gamma(shape,scale))
    void PriorResample(const Selector<double> &poswarray) {
        for (int i = 0; i < GetSize(); i++) {
            if (!poswarray.GetVal(i)) { (*this)[i] = Random::GammaSample(GetShape(), GetScale()); }
        }
    }

    //! get mean over the array
    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetSize(); i++) { m1 += GetVal(i); }
        m1 /= GetSize();
        return m1;
    }

    //! get variance over the array
    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetSize(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetSize();
        m2 /= GetSize();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const double &mean;
    const double &invshape;
};

template <>
struct has_custom_serialization<IIDGamma> : std::true_type {};

/**
 * \brief A tree-structured branch-wise array of iid Gamma variables
 * (tree-stuctured version of IIDGamma)
 */

class BranchIIDGamma : public SimpleBranchArray<double> {
  public:
    BranchIIDGamma(const Tree &intree, const double &inmean, const double &ininvshape)
        : SimpleBranchArray<double>(intree), mean(inmean), invshape(ininvshape) {
        Sample();
    }

    ~BranchIIDGamma() {}

    double GetShape() const { return 1.0 / invshape; }
    double GetScale() const { return 1.0 / invshape / mean; }

    //! set all entries equal to inval
    void SetAllBranches(double inval) {
        for (int i = 0; i < GetNbranch(); i++) { (*this)[i] = inval; }
    }

    //! sample all entries from prior
    void Sample() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] = Random::GammaSample(GetShape(), GetScale());
        }
    }

    //! resample all entries from posterior, conditional on BranchArray of
    //! PoissonSuffStat
    // void GibbsResample(const PoissonSuffStatBranchArray& suffstatarray)
    // {
    void GibbsResample(const BranchArray<PoissonSuffStat> &suffstatarray) {
        for (int i = 0; i < GetNbranch(); i++) {
            const PoissonSuffStat &suffstat = suffstatarray.GetVal(i);
            (*this)[i] = Random::GammaSample(
                GetShape() + suffstat.GetCount(), GetScale() + suffstat.GetBeta());
        }
    }

    //! get total log prob summed over all branches
    double GetLogProb() {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) { total += GetLogProb(i); }
        return total;
    }

    //! get log prob for a given branch
    double GetLogProb(int index) {
        return Random::logGammaDensity(GetVal(index), GetShape(), GetScale());
    }

    //! get sum over all entries (name is rather specialized... could change..)
    double GetTotalLength() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) { m1 += GetVal(i); }
        return m1;
    }

    //! get mean over the array
    double GetMean() const {
        double m1 = 0;
        for (int i = 0; i < GetNbranch(); i++) { m1 += GetVal(i); }
        m1 /= GetNbranch();
        return m1;
    }

    //! get variance over the array
    double GetVar() const {
        double m1 = 0;
        double m2 = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            m1 += GetVal(i);
            m2 += GetVal(i) * GetVal(i);
        }
        m1 /= GetNbranch();
        m2 /= GetNbranch();
        m2 -= m1 * m1;
        return m2;
    }

  protected:
    const double &mean;
    const double &invshape;
};

template <>
struct has_custom_serialization<BranchIIDGamma> : std::true_type {};
