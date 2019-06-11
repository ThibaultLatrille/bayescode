#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

class Tree;

/**
 * \brief A set of taxa (simply, a bi-directional correspondance between taxon
 * names and integers between 0 and Ntaxa-1)
 *
 * TaxonSet implements a globally consistent correspondance of taxon names and
 * indices between trees and sequence alignments. Each Tree and
 * SequenceAlignment object has a pointer to a TaxonSet. Typically, a
 * SequenceAlignment is created by reading from a file; the TaxonSet is created
 * at that step, by the SequenceAlignment constructor itself. Then, a Tree is
 * read from file and is registered with the taxon set given by the sequence
 * alignment. Once this is done, both objects have mutually consistent indexing
 * systems and can work together (e.g. to calculate a likelihood).
 *
 * Of note: some work to do here to ensure proper destruction of TaxonSet (smart
 * pointers?)
 */

class TaxonSet {
  public:
    //! constructor, based on a vector of taxon names
    explicit TaxonSet(const std::vector<std::string> &names);

    //! constructor, based on a vector of taxon names
    TaxonSet(const TaxonSet &from);

    //! default constructor
    ~TaxonSet() = default;

    //! return number of taxa
    int GetNtaxa() const { return Ntaxa; }

    //! return taxon name, given the index
    std::string GetTaxon(int index) const { return taxlist[index]; }

    //! return true if taxon is present, given the name
    bool TaxonPresent(std::string const &intaxon) const;

    //! return taxon index, given the name
    int GetTaxonIndex(std::string const &intaxon) const;

    //! return taxon index, given incomplete name (first part)
    int GetTaxonIndexWithIncompleteName(std::string const &taxname) const;

    //! formatted output to stream
    void ToStream(std::ostream &os) const;

    //! The taxon given the node
    std::vector<int> get_index_table(const Tree *tree) const;

    //! The node given the taxon
    std::vector<int> get_reverse_index_table(const Tree *tree) const;

  private:
    int Ntaxa;
    std::map<std::string, int> taxmap;
    std::vector<std::string> taxlist;
};

inline int TaxonSet::GetTaxonIndex(std::string const &intaxon) const {
    auto i = taxmap.find(intaxon);
    if (i == taxmap.end()) {
        std::cerr << "error in TaxonSet: taxon not found\n";
        exit(1);
    }
    return i->second - 1;
}

inline bool TaxonSet::TaxonPresent(std::string const &intaxon) const {
    return taxmap.find(intaxon) != taxmap.end();
}