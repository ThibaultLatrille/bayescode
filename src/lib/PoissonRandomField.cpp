#include "PoissonRandomField.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>

using namespace std;

PoissonRandomField::PoissonRandomField(
    set<unsigned> const &sample_size_set, CodonStateSpace const &instatespace, unsigned precision)
    : statespace{instatespace}, precision{precision} {
    unsigned nbr_steps = PowUnsigned(2, precision) + 1;
    std::cerr << "Number of steps (precision) of PRF: " << nbr_steps << std::endl;
    grid_s_step = 20.0 / nbr_steps;

    for (unsigned sample_size : sample_size_set) {
        ComputedProb[sample_size] = deque<pair<double, vector<double>>>();
        ComputedBinom[sample_size] = BinomialCoefficientArray(sample_size);

        ComputedProb.at(sample_size).push_back(make_pair(0, ExpectedTimeObsVector(sample_size, 0)));
    }
}

double PoissonRandomField::GetProb(int anc_state, int der_state, unsigned der_occurence,
    unsigned sample_size, const vector<double> &aafitnessarray, const GTRSubMatrix &nucmatrix,
    const double &theta) const {
    if (anc_state < 0 or der_state < 0) { return 0.0; }

    double proba_obs = 0;

    if (anc_state == der_state) {
        // If the ancestral allele is monomorphic

        proba_obs = 1.0;
        for (auto neighbor_state : statespace.GetNeighbors(anc_state)) {
            // 0 is the special case for which it is the sum of all over 0 < i <= n
            // Nucleotide mutation rate between ancestral and derived codon
            proba_obs -= theta * InterpolateProba(anc_state, neighbor_state, 0, sample_size,
                                     aafitnessarray, nucmatrix);
        }

    } else {
        // If the ancestral allele is not monomorphic
        proba_obs = theta * InterpolateProba(anc_state, der_state, der_occurence, sample_size,
                                aafitnessarray, nucmatrix);
    }
    assert(!std::isnan(proba_obs));
    if (proba_obs >= 0.0 and proba_obs <= 1.0) {
        return proba_obs;
    } else {
        return 0.0;
    };
}

static struct {
    bool operator()(
        const pair<double, vector<double>> &left, const pair<double, vector<double>> &right) {
        return left.first < right.first;
    }

    bool operator()(const pair<double, vector<double>> &left, float right) {
        return left.first < right;
    }

    bool operator()(float left, const pair<double, vector<double>> &right) {
        return left < right.first;
    }
} PairLowerThan;

double PoissonRandomField::InterpolateProba(int anc_state, int der_state, unsigned der_occurence,
    unsigned sample_size, const vector<double> &aafitnessarray, const GTRSubMatrix &nucmatrix) const {
    // Selection coefficient between ancestral and derived codon
    double s = log(aafitnessarray.at(statespace.Translation(der_state)));
    s -= log(aafitnessarray.at(statespace.Translation(anc_state)));

    if (ComputedProb.at(sample_size).front().first > s or
        ComputedProb.at(sample_size).back().first < s) {
        UpdateComputed(sample_size, s);
    }

    // The index for the closest (lower and upper) selection coefficient for which pre-computation
    // is available
    auto it_up = upper_bound(
        ComputedProb.at(sample_size).begin(), ComputedProb.at(sample_size).end(), s, PairLowerThan);
    auto it_low = prev(it_up);

    double f;
    if (it_low->first == s) {
        f = it_low->second.at(der_occurence);
    } else {
        // Linear interpolation using the closest (lower and upper) pre-computation available
        double p = (s - it_low->first) / (it_up->first - it_low->first);
        f = p * it_up->second.at(der_occurence) + (1 - p) * it_low->second.at(der_occurence);
    }


    int pos = statespace.GetDifferingPosition(anc_state, der_state);
    assert(0 <= pos and pos < 3);
    double mutation_rate = nucmatrix(
        statespace.GetCodonPosition(pos, anc_state), statespace.GetCodonPosition(pos, der_state));

    return mutation_rate * f;
}

void PoissonRandomField::UpdateComputed(unsigned sample_size, double s) const {
    while (ComputedProb.at(sample_size).front().first > s) {
        double new_s = ComputedProb.at(sample_size).front().first - grid_s_step;
        ComputedProb.at(sample_size)
            .push_front(make_pair(new_s, ExpectedTimeObsVector(sample_size, new_s)));
    }
    while (ComputedProb.at(sample_size).back().first < s) {
        double new_s = ComputedProb.at(sample_size).back().first + grid_s_step;
        ComputedProb.at(sample_size)
            .push_back(make_pair(new_s, ExpectedTimeObsVector(sample_size, new_s)));
    }
}

vector<double> PoissonRandomField::ExpectedTimeObsVector(unsigned n, double s) const {
    // n is the sample size
    // s is the selection coefficient associated to the derived allele

    vector<double> obs_array(n + 1, 0);

    // The grid on which we calculate the integral
    unsigned grid_size = PowUnsigned(2, precision);

    // The bounds of the integral
    double x_min = 0.0;
    double x_max = 1;

    // The step to increment x
    double h = (x_max - x_min) / grid_size;
    double x = x_min;

    // Intermediate calculus so they don't need to be recomputed n times
    vector<double> res_array(grid_size, 0);

    for (unsigned i{0}; i < grid_size; i++) {
        if (s > 1e-8) {
            res_array[i] = 2 * (1 - exp(-s * (1 - x))) / (1 - exp(-s));
        } else {
            res_array[i] = 2 * (1 - x);
        }
        x += h;
    }

    for (unsigned a{1}; a <= n; a++) {
        // Integral using trapezoidal rule
        double integral = 0;
        x = x_min;

        // x = 0 (trapezoidal rule so must be divided by 2)
        integral += h * res_array[0] * pow(x, a - 1) * pow(1 - x, n - a - 1) / 2;
        x += h;

        // 0 < x < 1
        for (unsigned i{1}; i < grid_size; i++) {
            integral += h * res_array[i] * pow(x, a - 1) * pow(1 - x, n - a - 1);
            x += h;
        }

        // x = 1, serie expansion at first order (trapezoidal rule so must be divided by 2)
        assert(abs(x - 1) < 1e-8);
        if (s > 1e-8) {
            integral += h * (s / (1 - exp(-s))) * pow(x, a - 1) * pow(1 - x, n - a);
        } else {
            integral += h * pow(x, a - 1) * pow(1 - x, n - a);
        }
        obs_array[a] = ComputedBinom.at(n).at(a) * integral;
    }

    // 0 is the special case for which it is the sum of all over 0 < i <= n
    assert(obs_array[0] == 0);
    obs_array[0] = SumVector(obs_array);
    return obs_array;
}


template <class T>
string join(vector<T> const &v, char sep) {
    return accumulate(v.begin() + 1, v.end(), to_string(v[0]),
        [sep](const string &acc, T b) { return acc + sep + to_string(b); });
};

vector<unsigned long long> BinomialCoefficientArray(unsigned n) {
    vector<unsigned long long> bin_array(n + 1, 1);

    for (unsigned k = 1; k <= n - k; k++) {
        bin_array[k] = bin_array[k - 1] / k;
        bin_array[k] *= n - k + 1;
        bin_array[n - k] = bin_array[k];
    }

    return bin_array;
}

unsigned PowUnsigned(unsigned a, unsigned b) {
    if (b == 0) return 1;
    if (b == 1) return a;

    unsigned tmp = PowUnsigned(a, b / 2);
    if (b % 2 == 0)
        return tmp * tmp;
    else
        return a * tmp * tmp;
}

double SumVector(vector<double> const &v) { return accumulate(v.begin(), v.end(), 0.0); }