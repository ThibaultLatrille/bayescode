#pragma once

#include <cassert>
#include <typeinfo>
#include "CodonSubMatrixArray.hpp"
#include "PathSuffStat.hpp"
#include "PoissonSuffStat.hpp"

/**
 * \brief A sufficient statistic for substitution histories, as a function of
 * the underlying nucleotide rate matrix (for Muse-Gaut codon models only)
 *
 * The generic sufficient statistics for substitution histories (S) as a
 * function of the rate matrix Q are defined in PathSuffStat. They give all
 * information needed in order to compute p(S | Q), up to a normalization
 * constant. When Q itself is a Muse and Gaut codon model parameterized by a
 * nucleotide rate matrix M, then the probability of S as a function of M, p(S |
 * M), can be expressed in terms of an even more compact (4x4) suff stat.
 *
 * NucPathSuffStat implements this idea, and provides methods for calculating
 * these suffstats based on generic PathSuffStat at the codon-level, adding them
 * over sites / branches and computing the log p(S | M) for any M.
 */

class NucPathSuffStat : public SuffStat {
  public:
    NucPathSuffStat()
        : rootcount(4, 0),
          paircount(4, std::vector<int>(4, 0)),
          pairbeta(4, std::vector<double>(4, 0)) {}

    ~NucPathSuffStat() {}

    //! set suff stat to 0
    void Clear() {
        for (int i = 0; i < Nnuc; i++) {
            rootcount[i] = 0;
            for (int j = 0; j < Nnuc; j++) {
                paircount[i][j] = 0;
                pairbeta[i][j] = 0;
            }
        }
    }

    //! \brief compute the 4x4 path suff stat out of 61x61 codonpathsuffstat
    //
    //! Note that the resulting 4x4 nuc path suff stat depends on other aspects of
    //! the codon matrix (e.g. the value of omega)
    void AddSuffStat(const NucCodonSubMatrix &codonmatrix, const PathSuffStat &codonpathsuffstat) {
        const CodonStateSpace *cod = codonmatrix.GetCodonStateSpace();
        const SubMatrix *nucmatrix = codonmatrix.GetNucMatrix();

        // root part
        const std::map<int, int> &codonrootcount = codonpathsuffstat.GetRootCountMap();
        for (std::map<int, int>::const_iterator i = codonrootcount.begin();
             i != codonrootcount.end(); i++) {
            int codon = i->first;
            rootcount[cod->GetCodonPosition(0, codon)] += i->second;
            rootcount[cod->GetCodonPosition(1, codon)] += i->second;
            rootcount[cod->GetCodonPosition(2, codon)] += i->second;
        }

        const std::map<int, double> &waitingtime = codonpathsuffstat.GetWaitingTimeMap();
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            int codon = i->first;
            for (int c2 = 0; c2 < cod->GetNstate(); c2++) {
                if (c2 != codon) {
                    int pos = cod->GetDifferingPosition(codon, c2);
                    if (pos < 3) {
                        int n1 = cod->GetCodonPosition(pos, codon);
                        int n2 = cod->GetCodonPosition(pos, c2);
                        pairbeta[n1][n2] +=
                            i->second * codonmatrix(codon, c2) / (*nucmatrix)(n1, n2);
                    }
                }
            }
        }

        const std::map<std::pair<int, int>, int> &codonpaircount =
            codonpathsuffstat.GetPairCountMap();
        for (std::map<std::pair<int, int>, int>::const_iterator i = codonpaircount.begin();
             i != codonpaircount.end(); i++) {
            int cod1 = i->first.first;
            int cod2 = i->first.second;
            int pos = cod->GetDifferingPosition(cod1, cod2);
            if (pos == 3) {
                std::cerr << "error in codon conj path suffstat\n";
                exit(1);
            }
            int n1 = cod->GetCodonPosition(pos, cod1);
            int n2 = cod->GetCodonPosition(pos, cod2);
            paircount[n1][n2] += i->second;
        }
    }

    //! \brief compute the 4x4 path suff stats out of 61x61 codonpathsuffstats
    //! over a one-dim array, and sum them up into this
    //!
    //! Note that this function assumes that each codonpathsuffstat given by the
    //! second array has a potentially different codon matrix, such as specified
    //! by the first array.
    void AddSuffStat(const MGOmegaCodonSubMatrixArray &codonmatrixarray,
        const PathSuffStatArray &codonpathsuffstatarray) {
        for (int i = 0; i < codonmatrixarray.GetSize(); i++) {
            AddSuffStat(codonmatrixarray.GetVal(i), codonpathsuffstatarray.GetVal(i));
        }
    }

    void AddSuffStat(const BranchSelector<MGOmegaCodonSubMatrix> &codonmatrixtree,
        const MGOmegaCodonSubMatrix &rootcodonmatrix,
        const PathSuffStatNodeArray &codonpathsuffstatnodearray) {
        // void AddSuffStat(const MGOmegaCodonSubMatrixBranchArray& codonmatrixtree,
        // const MGOmegaCodonSubMatrix& rootcodonmatrix, const
        // PathSuffStatNodeArray& codonpathsuffstatnodearray)    {
        RecursiveAddSuffStat(codonmatrixtree.GetTree(), codonmatrixtree.GetTree().root(),
            codonmatrixtree, rootcodonmatrix, codonpathsuffstatnodearray);
    }

    void RecursiveAddSuffStat(const Tree &tree, Tree::NodeIndex from,
        const BranchSelector<MGOmegaCodonSubMatrix> &codonmatrixtree,
        const MGOmegaCodonSubMatrix &rootcodonmatrix,
        const PathSuffStatNodeArray &codonpathsuffstatnodearray) {
        if (tree.is_root(from)) {
            AddSuffStat(rootcodonmatrix, codonpathsuffstatnodearray.GetVal(from));
        } else {
            AddSuffStat(codonmatrixtree.GetVal(tree.branch_index(from)),
                codonpathsuffstatnodearray.GetVal(from));
        }
        for (auto c : tree.children(from)) {
            RecursiveAddSuffStat(
                tree, c, codonmatrixtree, rootcodonmatrix, codonpathsuffstatnodearray);
        }
    }

    //! \brief return the log probability as a function of a nucleotide matrix
    //!
    //! The codon state space is given as an argument (the nucleotide matrix or
    //! the suff stat themselves do not know the genetic code)
    double GetLogProb(const SubMatrix &mat, const CodonStateSpace &cod) const {
        double total = 0;
        // root part
        int nroot = 0;
        auto rootstat = mat.GetStationary();
        for (int i = 0; i < Nnuc; i++) {
            total += rootcount[i] * log(rootstat[i]);
            nroot += rootcount[i];
        }
        total -= nroot / 3 * log(cod.GetNormStat(rootstat));

        // non root part
        for (int i = 0; i < Nnuc; i++) {
            for (int j = 0; j < Nnuc; j++) {
                if (i != j) {
                    total += paircount[i][j] * log(mat(i, j));
                    total -= pairbeta[i][j] * mat(i, j);
                }
            }
        }

        return total;
    }

    //! add another nucpath suff stat to this
    void Add(const NucPathSuffStat &from) {
        for (int i = 0; i < Nnuc; i++) { rootcount[i] += from.rootcount[i]; }
        for (int i = 0; i < Nnuc; i++) {
            for (int j = 0; j < Nnuc; j++) {
                paircount[i][j] += from.paircount[i][j];
                pairbeta[i][j] += from.pairbeta[i][j];
            }
        }
    }

    //! add another nuc path suffstat to this, operator version
    NucPathSuffStat &operator+=(const NucPathSuffStat &from) {
        Add(from);
        return *this;
    }

    template <class T>
    void serialization_interface(T &x) {
        x.add(rootcount, paircount, pairbeta);
    }

    // private:
    std::vector<int> rootcount;
    std::vector<std::vector<int>> paircount;
    std::vector<std::vector<double>> pairbeta;
};

template <>
struct has_custom_serialization<NucPathSuffStat> : std::true_type {};

/**
 * \brief A sufficient statistic for substitution histories, as a function of
 * omega=dN/dS (for codon models)
 *
 * The generic sufficient statistics for substitution histories (S) as a
 * function of the rate matrix Q are defined in PathSuffStat. They give all
 * information needed in order to compute p(S | Q), up to a normalization
 * constant. When Q itself is codon model with an omega=dN/dS acting as a
 * multiplier in front of all non-synonymous substitutions, then the probability
 * of S as a function of omega can be expressed in very compact form: p(S |
 * omega) propto omega^count exp(-beta * omega), where count (integer) and beta
 * (positive real number) are the suff stats. Note that this is in fact
 * analogous to a Poisson distribution, with mean omega, and thus,
 * OmegaPathSuffStat derives from PoissonSuffStat.
 */

class OmegaPathSuffStat : public PoissonSuffStat {
  public:
    OmegaPathSuffStat() {}
    ~OmegaPathSuffStat() {}

    //! \brief tease out syn and non-syn substitutions and sum up count and beta
    //! stats from a 61x61 codon path suffstat
    //!
    //! note that omega suff stat depends on the other aspects of the codon matrix
    //! (in particular, the nucleotide rate matrix)
    void AddSuffStat(const OmegaCodonSubMatrix &codonsubmatrix, const PathSuffStat &pathsuffstat) {
        int ncodon = codonsubmatrix.GetNstate();
        const CodonStateSpace *statespace = codonsubmatrix.GetCodonStateSpace();

        const std::map<std::pair<int, int>, int> &paircount = pathsuffstat.GetPairCountMap();
        const std::map<int, double> &waitingtime = pathsuffstat.GetWaitingTimeMap();

        double tmpbeta = 0;
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            double totnonsynrate = 0;
            int a = i->first;
            for (int b = 0; b < ncodon; b++) {
                if (b != a) {
                    if (codonsubmatrix(a, b) != 0) {
                        if (!statespace->Synonymous(a, b)) {
                            totnonsynrate += codonsubmatrix(a, b);
                        }
                    }
                }
            }
            tmpbeta += i->second * totnonsynrate;
        }
        tmpbeta /= codonsubmatrix.GetOmega();

        int tmpcount = 0;
        for (std::map<std::pair<int, int>, int>::const_iterator i = paircount.begin();
             i != paircount.end(); i++) {
            if (!statespace->Synonymous(i->first.first, i->first.second)) { tmpcount += i->second; }
        }

        PoissonSuffStat::AddSuffStat(tmpcount, tmpbeta);
    }

    //! \brief get count and beta stats from an array of 61x61 codon path suffstat
    //! and sum them up into this suff stat
    //!
    //! This method assumes that each codonpathsuffstat given by the second array
    //! has a potentially different codon matrix, such as specified by the first
    //! array.
    void AddSuffStat(const Selector<AAMutSelOmegaCodonSubMatrix> &codonsubmatrixarray,
        const Selector<PathSuffStat> &pathsuffstatarray) {
        for (int i = 0; i < codonsubmatrixarray.GetSize(); i++) {
            AddSuffStat(codonsubmatrixarray.GetVal(i), pathsuffstatarray.GetVal(i));
        }
    }

    double GetLogProb(double omega) const { return count * log(omega) - beta * omega; }
};

/**
 * \brief An array of omega suff stats
 */

class OmegaPathSuffStatArray : public SimpleArray<OmegaPathSuffStat>,
                               public Array<PoissonSuffStat> {
  public:
    //! constructor (param: array size)
    OmegaPathSuffStatArray(int insize) : SimpleArray<OmegaPathSuffStat>(insize) {}

    //! need to redefine GetSize(), because of mutiple inheritance, as an array of
    //! PoissonSuffstat and OmegaPathSuffStat
    int GetSize() const /*override*/ { return array.size(); }

    //! need to re-define GetVal(), by explicitly returning a const
    //! OmegaPathSuffStat&, because of multiple inheritance
    const OmegaPathSuffStat &GetVal(int i) const /*override*/ { return array[i]; }

    //! need to re-define operator[], by explicitly returning an
    //! OmegaPathSuffStat&, because of multiple inheritance
    OmegaPathSuffStat &operator[](int i) /*override*/ { return array[i]; }

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Clear(); }
    }

    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    void AddSuffStat(const Selector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
        const Selector<PathSuffStat> &pathsuffstatarray) {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].AddSuffStat(codonsubmatrixarray.GetVal(i), pathsuffstatarray.GetVal(i));
        }
    }

    void AddSuffStat(const Selector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
        const NodeSelector<PathSuffStat> &pathsuffstatnodearray, const BranchSelector<int> &alloc) {
        RecursiveAddSuffStat(alloc.GetTree(), alloc.GetTree().root(), codonsubmatrixarray,
            pathsuffstatnodearray, alloc);
    }

    void RecursiveAddSuffStat(const Tree &tree, Tree::NodeIndex from,
        const Selector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
        const NodeSelector<PathSuffStat> &pathsuffstatnodearray, const BranchSelector<int> &alloc) {
        if (!tree.is_root(from)) {
            int i = alloc.GetVal(from);
            (*this)[i].AddSuffStat(
                codonsubmatrixarray.GetVal(i), pathsuffstatnodearray.GetVal(from));
        }
        for (auto c : tree.children(from)) {
            RecursiveAddSuffStat(tree, c, codonsubmatrixarray, pathsuffstatnodearray, alloc);
        }
    }

    //! compute omega suff stats and do a member-wise addition -- for
    //! mutation-selection codon matrices
    void AddSuffStat(const Selector<AAMutSelOmegaCodonSubMatrix> &codonsubmatrixarray,
        const Selector<PathSuffStat> &pathsuffstatarray) {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].AddSuffStat(codonsubmatrixarray.GetVal(i), pathsuffstatarray.GetVal(i));
        }
    }

    //! \brief add suffstatarray given as argument to this array based on the
    //! allocations provided as the second argument (mixture models)
    //!
    //! specifically, for each i=0..GetSize()-1, (*this)[alloc[i]] +=
    //! suffstatarray[i]
    void Add(const Selector<OmegaPathSuffStat> &suffstatarray, const Selector<int> &alloc) {
        assert(suffstatarray.GetSize() == alloc.GetSize());
        for (int i = 0; i < suffstatarray.GetSize(); i++) {
            (*this)[alloc.GetVal(i)].Add(suffstatarray.GetVal(i));
        }
    }

    //! return total log prob over array, given an array of omega_i's of same size
    double GetLogProb(const Array<double> *omegaarray) const {
        assert(omegaarray->GetSize() == GetSize());
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetVal(i).GetLogProb(omegaarray->GetVal(i));
        }
        return total;
    }

    //! return marginal log prob summed over omega_i's that are iid gamma(shape,
    //! scale) (conjugate gamma-poisson relation)
    double GetMarginalLogProb(double shape, double scale) const {
        double total = 0;
        // factoring out prior factor
        for (int i = 0; i < GetSize(); i++) {
            int count = GetVal(i).GetCount();
            double beta = GetVal(i).GetBeta();
            total += -(shape + count) * log(scale + beta) + Random::logGamma(shape + count);
        }
        total += GetSize() * (shape * log(scale) - Random::logGamma(shape));
        return total;
    }
};

/**
 * \brief A BranchArray of omega suff stats
 */

class OmegaPathSuffStatBranchArray : public SimpleBranchArray<OmegaPathSuffStat>,
                                     public BranchArray<PoissonSuffStat> {
  public:
    //! constructor (param: tree)
    OmegaPathSuffStatBranchArray(const Tree &intree)
        : SimpleBranchArray<OmegaPathSuffStat>(intree), tree(intree) {}
    ~OmegaPathSuffStatBranchArray() {}

    //! need to redefine GetTree(), because of mutiple inheritance, as an array of
    //! PoissonSuffstat and OmegaPathSuffStat
    const Tree &GetTree() const /*override*/ { return tree; }

    int GetNbranch() const { return GetTree().nb_branches(); }

    //! need to re-define GetVal(), by explicitly returning a const
    //! OmegaPathSuffStat&, because of multiple inheritance
    const OmegaPathSuffStat &GetVal(int i) const /*override*/ { return array[i]; }

    //! need to re-define operator[], by explicitly returning an
    //! OmegaPathSuffStat&, because of multiple inheritance
    OmegaPathSuffStat &operator[](int i) /*override*/ { return array[i]; }

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetNbranch(); i++) { (*this)[i].Clear(); }
    }

    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    void AddSuffStat(const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
        const MGOmegaCodonSubMatrix &rootcodonsubmatrix,
        const NodeSelector<PathSuffStat> &pathsuffstatarray) {
        RecursiveAddSuffStat(
            GetTree().root(), codonsubmatrixarray, rootcodonsubmatrix, pathsuffstatarray);
    }

    void RecursiveAddSuffStat(Tree::NodeIndex from,
        const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
        const MGOmegaCodonSubMatrix &rootcodonsubmatrix,
        const NodeSelector<PathSuffStat> &pathsuffstatarray) {
        if (!GetTree().is_root(from)) {
            (*this)[GetTree().branch_index(from)].AddSuffStat(
                codonsubmatrixarray.GetVal(GetTree().branch_index(from)),
                pathsuffstatarray.GetVal(from));
        }
        for (auto c : GetTree().children(from)) {
            RecursiveAddSuffStat(c, codonsubmatrixarray, rootcodonsubmatrix, pathsuffstatarray);
        }
    }

    //! return total log prob over array, given an array of omega_i's of same size
    double GetLogProb(const BranchArray<double> *omegaarray) const {
        double total = 0;
        for (int i = 0; i < GetNbranch(); i++) {
            total += GetVal(i).GetLogProb(omegaarray->GetVal(i));
        }
        return total;
    }

    //! return marginal log prob summed over omega_i's that are iid gamma(shape,
    //! scale) (conjugate gamma-poisson relation)
    double GetMarginalLogProb(double shape, double scale) const {
        double total = 0;
        // factoring out prior factor
        for (int i = 0; i < GetNbranch(); i++) {
            int count = GetVal(i).GetCount();
            double beta = GetVal(i).GetBeta();
            total += -(shape + count) * log(scale + beta) + Random::logGamma(shape + count);
        }
        total += GetNbranch() * (shape * log(scale) - Random::logGamma(shape));
        return total;
    }

  private:
    const Tree &tree;
};
