#pragma once

#include <cassert>
#include "components/Tracer.hpp"

template <class Model>
bool check_restart(Model& model, std::string tracefile) {
    // get stats from new model
    Tracer tracer(model, processing::HasTag<Stat>());
    std::stringstream ss;
    tracer.write_line(ss);

    // get stats from trace file
    std::ifstream trace{tracefile};
    std::string line, nonempty_line;
    while (std::getline(trace, line)) {
        if (line != "\n") { nonempty_line = line; }
    }

    if (("\n" + nonempty_line) != ss.str()) {
        std::cerr << "Warning: statistics from restarted model do not match statistics from "
                     "tracefile "
                  << tracefile << "\n";

        std::stringstream header;
        tracer.write_header(header);
        std::stringstream restarted;
        tracer.write_line(restarted);
        std::stringstream loaded(nonempty_line);

        std::string each;
        std::vector<std::string> header_vector;
        std::vector<std::string> restarted_vector;
        std::vector<std::string> loaded_vector;
        while (std::getline(header, each, '\t')) { header_vector.push_back(each); };
        while (std::getline(restarted, each, '\t')) { restarted_vector.push_back(each); };
        while (std::getline(loaded, each, '\t')) { loaded_vector.push_back(each); };

        assert(header_vector.size() == restarted_vector.size());
        assert(loaded_vector.size() == restarted_vector.size());

        // Remove "\n" at the beginning of restarted trace
        restarted_vector[0].erase(0, 1);

        std::cout << "Statistic\tLoaded\tRestarded\t" << std::endl;
        for (size_t word = 0; word < header_vector.size(); word++) {
            bool diff = loaded_vector[word] != restarted_vector[word];
            if (diff) { std::cout << "\e[35m\e[1m"; }
            std::cout << header_vector[word] << "\t" << loaded_vector[word] << "\t"
                      << restarted_vector[word] << "\t" << std::endl;
            if (diff) { std::cout << "\e[0m"; }
        }
        return false;
    } else {
        std::cout << "Restart check passed successfully!\n";
        std::cout << "Statistics are:\n";
        tracer.write_header(std::cout);
        tracer.write_line(std::cout);
        std::cerr << "\n";
        return true;
    }
}