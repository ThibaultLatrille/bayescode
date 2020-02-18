#pragma once

#include <numeric>
#include "components/ChainReader.hpp"
#include "components/Tracer.hpp"

// Mean of a vector
double mean(std::vector<double> const &v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

// Variance of a vector
double var(std::vector<double> const &v) {
    double s2 = std::accumulate(v.begin(), v.end(), 0.0, [](double a, double const &b) {
        return a + b * b;
    }) / v.size();
    double s = mean(v);
    return s2 - s * s;
}

template <class Model>
void stats_posterior(Model &model, ChainReader &cr, int const &every, int const &size) {
    unsigned nbr_header_fields{0};
    std::stringstream ss_header;
    {
        Tracer tracer(model, processing::HasTag<Stat>());
        tracer.write_header(ss_header);
        nbr_header_fields = tracer.nbr_header_fields();
    }

    std::vector<std::vector<double>> bidim_stats(nbr_header_fields, std::vector<double>(size, 0));

    for (int step = 0; step < size; step++) {
        std::cerr << '.';
        cr.skip(every);
        Tracer tracer(model, processing::HasTag<Stat>());
        model.Update();
        std::vector<double> line_values = tracer.line_values();
        for (unsigned field{0}; field < nbr_header_fields; field++) {
            bidim_stats.at(field).at(step) = line_values.at(field);
        }
    }
    std::cerr << '\n';

    std::string stat;
    unsigned field{0};
    while (getline(ss_header, stat, '\t')) {
        std::cout << "posterior " << stat << " : " << mean(bidim_stats.at(field)) << "±"
                  << sqrt(var(bidim_stats.at(field))) << '\n';
        field++;
    }
}

template <class Model>
void recompute_trace(
    Model &model, ChainReader &cr, std::string const &name, int const &every, int const &size) {
    Tracer tracer(model, processing::HasTag<Stat>());
    std::ofstream os(name + ".trace.tsv");
    tracer.write_header(os);

    for (int step = 0; step < size; step++) {
        std::cerr << '.';
        cr.skip(every);
        model.Update();
        tracer.write_line(os);
    }
    std::cerr << '\n';
}