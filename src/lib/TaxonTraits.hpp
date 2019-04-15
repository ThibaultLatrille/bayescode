#pragma once

#include <fstream>
#include <sstream>
#include "Random.hpp"
#include "TaxonMapping.hpp"

class TaxonTraits {
  public:
    //! constructor, based on a vector of taxon names
    explicit TaxonTraits(
        std::string const &filename, TaxonSet const &taxon_set, bool polymorphism_aware) {
        taxon_presence.resize(taxon_set.GetNtaxa(), false);
        int gentime_dim{-1};

        // for line in file
        std::ifstream input_stream{filename};
        if (!input_stream) {
            std::cerr << "Traits file " << filename << " doesn't exist" << std::endl;
        }

        std::string line;
        // skip the header of the file
        getline(input_stream, line);
        char sep{'\t'};
        {
            std::istringstream line_stream(line);
            std::string word;
            // Skip the first column (taxon name)
            getline(line_stream, word, sep);
            assert(word == "TaxonName");
            int counter = 0;
            while (getline(line_stream, word, sep)) {
                if (word == "LogGenerationTime" and polymorphism_aware) {
                    gentime = true;
                    taxon_gentime_presence.resize(taxon_set.GetNtaxa(), false);
                    taxon_gentime.resize(taxon_set.GetNtaxa(), 0.0);
                    gentime_dim = counter;
                } else {
                    header.push_back(word);
                }
                counter++;
            }
        }
        dimensions = header.size();
        if (!header.empty()) {
            taxon_traits_presence.resize(taxon_set.GetNtaxa());
            taxon_traits.resize(taxon_set.GetNtaxa());
        }
        while (getline(input_stream, line)) {
            std::istringstream line_stream(line);
            std::string taxon{};
            getline(line_stream, taxon, sep);
            int id = taxon_set.GetTaxonIndex(taxon);
            taxon_presence.at(id) = true;
            if (!header.empty()) {
                taxon_traits_presence[id] = std::vector<bool>(dimensions, false);
                taxon_traits[id] = std::vector<double>(dimensions, 0.0);
            }
            std::string word;
            int counter = 0;
            int dim_counter = 0;
            while (getline(line_stream, word, sep)) {
                if (!word.empty() and convertible_to_float(word)) {
                    if (counter == gentime_dim) {
                        taxon_gentime_presence.at(id) = true;
                        taxon_gentime.at(id) = std::stod(word);
                    } else {
                        taxon_traits_presence.at(id).at(dim_counter) = true;
                        taxon_traits.at(id).at(dim_counter) = std::stod(word);
                    }
                }
                if (counter != gentime_dim) { dim_counter++; }
                counter++;
            }
        }
        assert(polymorphism_aware == gentime);
    }

    bool static convertible_to_float(const string &in) {
        std::stringstream sstr(in);
        float f;
        return bool(sstr >> f);
    }

    int GetDim() const { return dimensions; }
    int TraitDimToMultivariateDim(int dim) const { return 2 + dim + gentime; }

    bool DataPresence(int taxon, int dim) const {
        if (taxon_presence[taxon]) {
            return taxon_traits_presence[taxon][dim];
        } else {
            return false;
        }
    }

    double Data(int taxon, int dim) const {
        assert(taxon_traits_presence[taxon][dim]);
        return taxon_traits[taxon][dim];
    }

    bool GenTimePresence() const { return gentime; }

    bool GenTimePresence(int taxon) const { return taxon_gentime_presence[taxon]; }

    double GenTime(int taxon) const {
        assert(taxon_gentime_presence[taxon]);
        return taxon_gentime[taxon];
    }

    //! default constructor
    ~TaxonTraits() = default;

  private:
    int dimensions{0};
    std::vector<std::string> header;
    std::vector<bool> taxon_presence;
    std::vector<std::vector<bool>> taxon_traits_presence;
    std::vector<std::vector<double>> taxon_traits;

    bool gentime{false};
    std::vector<bool> taxon_gentime_presence;
    std::vector<double> taxon_gentime;
};
