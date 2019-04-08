#pragma once

#include <map>
#include "Random.hpp"
#include "StateSpace.hpp"

/**
 * \brief A codon state space
 *
 * CodonStateSpace implements the methods of the base class StateSpace
 * (conversion from int to string and conversely -- here, the strings are
 * nucleotide triplets, e.g. "AAC"). By default, the list of codons considered
 * by any method or algorithm excludes stop codons. Thus, for instance, in the
 * case of the universal genetic code, Nstate==61. If a method takes or returns
 * codons including stops, then, this is made explicit in the method's name
 */

class CodonStateSpace : public StateSpace {
  public:
    static const int Npos = 3;

    //! constructor: should always specify the genetic code (en enum type:
    //! Universal, MtMam or MtInv, see BiologicalSequences.h)
    explicit CodonStateSpace(GeneticCodeType type);
    ~CodonStateSpace() throw() override;

    int GetNstate() const override { return Nstate; }

    //! return length of symbol used when printing state (normally, 1 for
    //! nucleotides or amino-acids, 3 for codons)
    virtual int GetSymbolLength() const override { return 3; }

    //! given a 3-nucleotide string, returns codon index, in 0..Nstate-1 (if stop
    //! or unrecognized, exits with error message)
    int GetState(std::string word) const override;

    //! given a codon index (stops excluded), returns the corresponding
    //! 3-nucleotide string
    std::string GetState(int codon) const override;

    // codon specific methods

    //! return the underlying nucleotide state space
    const DNAStateSpace *GetDNAStateSpace() const { return nucstatespace; }

    //! return the amino-acid state space (for translated sequences)
    const ProteinStateSpace *GetProteinStateSpace() const { return protstatespace; }

    //! \brief return a codon based on three nucleotides encoded as integers (see
    //! DNAStateSpace)
    //!
    //! returns -1 (== unknown) if at least one of the positions is unknown;
    //! if stop, then exits with error message.
    int GetCodonFromDNA(int pos1, int pos2, int pos3) const;

    //! \brief given 2 nearest-neighbor codons, returns at which position they
    //! differ
    //
    //! codons should not be stop codons;
    //! returns -1 if codons are identical;
    //! returns 3 if codons differ at more than one position;
    //! otherwise, returns the position at which codons differ (i.e. returns 0,1
    //! or 2 if the codons differ at position 1,2 or 3).
    int GetDifferingPosition(int i, int j) const;

    //! \brief return the vector of codons differing at exactly one position
    std::vector<int> GetNeighbors(int i) const;

    //! return the integer encoding for the nucleotide at requested position
    //! pos=0,1, or 2
    int GetCodonPosition(int pos, int codon) const {
        if ((pos < 0) || (pos >= Npos)) {
            std::cerr << "GetCodonPosition: pos out of bound\n";
            std::cerr << pos << '\n';
            exit(1);
        }
        if (codon == -1) { return -1; }
        if ((codon < 0) || (codon >= Nstate)) {
            std::cerr << "GetCodonPosition: codon out of bound\n";
            std::cerr << codon << '\n';
            exit(1);
        }
        return CodonPos[pos][codon];
    }

    //! translation: amino-acid encoded by given codon (stops excluded)
    int Translation(int codon) const { return CodonCode[codon]; }

    //! whether the two codons are synonymous or not
    bool Synonymous(int codon1, int codon2) const {
        return (CodonCode[codon1] == CodonCode[codon2]);
    }

    //! check whether the combination of the three integer-encoded nucleotides
    //! make a stop codon or not
    bool CheckStop(int pos1, int pos2, int pos3) const;

    //! computes the sum of nuc stats over stop codons (S) and returns 1-S
    double GetNormStat(const EVector &nucstat) const {
        double stopstat = 0;
        for (int i = 0; i < Nstop; i++) {
            stopstat += nucstat[StopPos1[i]] * nucstat[StopPos2[i]] * nucstat[StopPos3[i]];
        }
        return 1.0 - stopstat;
    }

    //! computes the sum of nuc stats over stop codons (S) and returns 1-S
    double GetNormStat(const double *nucstat) const {
        double stopstat = 0;
        for (int i = 0; i < Nstop; i++) {
            stopstat += nucstat[StopPos1[i]] * nucstat[StopPos2[i]] * nucstat[StopPos3[i]];
        }
        return 1.0 - stopstat;
    }

  private:
    // number of stop codons under this genetic code (typically 3 for the
    // Universal code)
    int GetNstop() const { return Nstop; }

    const int *GetStopPos1() const { return StopPos1; }

    const int *GetStopPos2() const { return StopPos2; }

    const int *GetStopPos3() const { return StopPos3; }

    GeneticCodeType code;
    const DNAStateSpace *nucstatespace;
    const ProteinStateSpace *protstatespace;
    // number of codons, not including stops (61 in general)
    int Nstate;

    // and array of size Ncodon = 64
    // whose entries are between -1 and 19
    // -1 : stop codon
    // 0..19 : amino acid encoded (1 letter code, alphabetical order)
    int *CodonCodeWithStops;
    int *CodonCode;
    int **CodonPos;
    int *StopCodons;
    int Nstop;
    int *StopPos1;
    int *StopPos2;
    int *StopPos3;

    mutable std::map<int, int> degeneracy;
    std::vector<std::vector<int>> neighbors_vector;
};