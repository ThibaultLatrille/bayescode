#include "PolySuffStat.hpp"

double PolySuffStat::GetLogProb(PoissonRandomField &poissonrandomfield,
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

void PolySuffStat::Add(const PolySuffStat &suffstat) {
    for (auto const &i : suffstat.GetPairCountMap()) { AddPairCount(i.first, i.second); }
}

int PolySuffStat::GetPairCount(std::tuple<int, int, unsigned, unsigned> poly_tuple) const {
    auto const i = polycount.find(poly_tuple);
    if (i == polycount.end()) { return 0; }
    return i->second;
}


void PolySuffStatArray::Clear() {
    for (int i = 0; i < GetSize(); i++) { (*this)[i].Clear(); }
}

double PolySuffStatArray::GetLogProb(PoissonRandomField &poissonrandomfield,
    const Selector<std::vector<double>> &siteaafitnessarray, const GTRSubMatrix &nucmatrix,
    const double &theta) const {
    double total = 0;
    for (int i = 0; i < GetSize(); i++) {
        total += GetVal(i).GetLogProb(
            poissonrandomfield, siteaafitnessarray.GetVal(i), nucmatrix, theta);
    }
    return total;
}

void PolySuffStatArray::Add(
    const Selector<PolySuffStat> &suffstatarray, const Selector<int> &alloc) {
    assert(alloc.GetSize() == suffstatarray.GetSize());
    for (int i = 0; i < suffstatarray.GetSize(); i++) {
        (*this)[alloc.GetVal(i)] += suffstatarray.GetVal(i);
    }
}

void PolySuffStatBidimArray::Clear() {
    for (int row = 0; row < GetNrow(); row++) {
        for (int col = 0; col < GetNcol(); col++) { (*this)(row, col).Clear(); }
    }
}

double PolySuffStatBidimArray::GetLogProb(PoissonRandomField &poissonrandomfield,
    const Selector<std::vector<double>> &siteaafitnessarray, const GTRSubMatrix &nucmatrix,
    const ScaledMutationRate &theta) const {
    double total = 0;
    for (int row = 0; row < GetNrow(); row++) {
        double d_theta = theta.GetTheta(row);
        for (int col = 0; col < GetNcol(); col++) {
            total += GetVal(row, col).GetLogProb(
                poissonrandomfield, siteaafitnessarray.GetVal(col), nucmatrix, d_theta);
        };
    }
    return total;
}

void PolySuffStatBidimArray::Add(
    const BidimSelector<PolySuffStat> &suffstatbidimarray, const Selector<int> &col_alloc) {
    assert(suffstatbidimarray.GetNrow() == this->GetNrow());
    assert(suffstatbidimarray.GetNcol() == col_alloc.GetSize());
    for (int row = 0; row < suffstatbidimarray.GetNrow(); row++) {
        for (int col = 0; col < suffstatbidimarray.GetNcol(); col++) {
            (*this)(row, col_alloc.GetVal(col)) += suffstatbidimarray.GetVal(row, col);
        }
    }
}
