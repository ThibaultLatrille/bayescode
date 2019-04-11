#pragma once

#include <vector>
#include "lib/SequenceAlignment.hpp"
#include "tree/interface.hpp"

class TaxonMap {
  public:
    //! constructor, based on a vector of taxon names
    TaxonMap(const Tree *tree, const SequenceAlignment *data) {
        taxon_table = data->GetTaxonSet()->get_index_table(tree);
        reverse_taxon_table = data->GetTaxonSet()->get_reverse_index_table(tree);
    };

    //! default constructor
    ~TaxonMap() = default;

    Tree::NodeIndex TaxonToNode(int taxon) const { return reverse_taxon_table[taxon]; }

    int NodeToTaxon(Tree::NodeIndex node) const {
        assert(taxon_table.at(node) != -1);
        return taxon_table[node];
    }

    int GetNnodes() const {
        return taxon_table.size();
    }

    int GetNtaxa() const {
        return reverse_taxon_table.size();
    }

  private:
    std::vector<int> taxon_table;
    std::vector<int> reverse_taxon_table;
};