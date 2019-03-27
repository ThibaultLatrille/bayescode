#pragma once

#include <string>
#include "BiologicalSequences.hpp"  //FIXME only here because constant unknown

/**
 * \brief Generic interface for a state space (nucleotide, codons, etc)
 *
 * In practice, states are encoded as integers ranging from 0 to Nstate-1, where
 * Nstate=4 for nucleotides, 61 for codons (universal code), etc. This is true
 * in particular for matrix calculation (SubMatrix), for sequence alignments
 * (SequenceAlignment). However, the StateSpace class provides the general
 * interface for converting integer states into strings or conversely.
 * Specialized versions, such as DNAStateSpace or CodonStateSpace provide
 * context-specific additional methods.
 */

class StateSpace {
  public:
    virtual ~StateSpace() = default;

    //! return size of state space
    virtual int GetNstate() const = 0;

    //! return integer for a given string-formatted state
    virtual int GetState(std::string from) const = 0;

    //! return a formatted string output for given integer-encoded state
    virtual std::string GetState(int state) const = 0;

    //! whether the two states are compatible (for the moment: just whether or not
    //! the two states are equal)
    virtual bool isCompatible(int state1, int state2) const {
        return ((state1 == unknown) || (state2 == unknown) || (state1 == state2));
    }

    //! return length of symbol used when printing state (normally, 1 for
    //! nucleotides or amino-acids, 3 for codons)
    virtual int GetSymbolLength() const { return 1; }
};

/**
 * \brief A StateSpace class built from an array of characters
 *
 * This class provides the implementation for all cases where states are
 * canonically referred to using a one-letter code (thus, nucleotides,
 * amino-acids, RNA). This excludes the case of codons.
 */

class SimpleStateSpace : public StateSpace {
  public:
    int GetState(std::string from) const override;

    int GetNstate() const override { return Nstate; }

    std::string GetState(int state) const override;

  protected:
    int Nstate;
    char *Alphabet;
    int NAlphabetSet;
    char *AlphabetSet;
};

/**
 * \brief DNA state space: A=0, C=1, G=2, T=3; Nstate=4.
 */

class DNAStateSpace : public SimpleStateSpace {
  public:
    DNAStateSpace();
    ~DNAStateSpace() throw() override;
};

/**
 * \brief RNA state space: A=0, C=1, G=2, U=3; Nstate=4.
 */

class RNAStateSpace : public SimpleStateSpace {
  public:
    RNAStateSpace();
    ~RNAStateSpace() throw() override;
};

/**
 * \brief Amino-acid state space
 *
 * Amino-acids are in alphabetical order based on the one-letter code:
 * ADEFGHIKLMNPQRSTVWY, Nstate=20. Note that this is not necessarily standard
 * (e.g. baseml or mrbayes use the alphabetical order for the three-letter
 * code).
 */

class ProteinStateSpace : public SimpleStateSpace {
  public:
    ProteinStateSpace();
    ~ProteinStateSpace() throw() override;
};

/**
 * \brief State space for RY (purine / pyrimidine) recoded nucleotide data
 */

class RYStateSpace : public SimpleStateSpace {
  public:
    RYStateSpace();
    ~RYStateSpace() throw() override;

    int GetRYCoding(int from);
};

/**
 * \brief A generic state space, for abitrary sets of characters (one-letter
 * code only)
 */

class GenericStateSpace : public SimpleStateSpace {
  public:
    GenericStateSpace(int inNstate, char *inAlphabet, int inNAlphabetSet, char *inAlphabetSet) {
        Nstate = inNstate;
        Alphabet = new char[Nstate];
        for (int i = 0; i < Nstate; i++) { Alphabet[i] = inAlphabet[i]; }
        NAlphabetSet = inNAlphabetSet;
        AlphabetSet = new char[NAlphabetSet];
        for (int i = 0; i < NAlphabetSet; i++) { AlphabetSet[i] = inAlphabetSet[i]; }
    }

    ~GenericStateSpace() throw() override {
        delete[] Alphabet;
        delete[] AlphabetSet;
    }
};