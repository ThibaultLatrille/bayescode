#pragma once

#include <fstream>
#include <map>
#include "BidimArray.hpp"
#include "BranchSitePath.hpp"
#include "BranchSiteSelector.hpp"
#include "Chrono.hpp"
#include "NodeArray.hpp"
#include "PolyProcess.hpp"
#include "SequenceAlignment.hpp"
#include "SubMatrix.hpp"
#include "tree/implem.hpp"

// PhyloProcess is a dispatcher:
// its responsibility is to create a random branch/site path
// for each branch/site pair

/**
 * \brief The core class of phylogenetic likelihood calculation and stochastic
 * mapping of substitution histories
 *
 * PhyloProcess takes as an input a tree, a sequence alignment (data), a set of
 * branch lengths, of site-specific rates and a selector of substitution
 * matrices across sites and branches. It is then responsible for organizing all
 * likelihood calculations by pruning, as stochastic mapping of substitution
 * histories. If polymorphism data is available (polyprocess is not a null pointer),
 * the likelihood of the data (number of occurrences in the population of the reference
 * and derived alleles at each site) is calculated using diffusion equations.
 */

class PhyloProcess {
  public:
    friend class PathSuffStat;
    friend class PathSuffStatArray;
    friend class PathSuffStatBidimArray;
    friend class PolySuffStat;
    friend class PolySuffStatArray;
    friend class PoissonSuffStatBranchArray;
    friend class PoissonSuffStatArray;
    friend class PathSuffStatNodeArray;

    //! \brief delegated constructor
    //!
    //! Constructor takes as parameters (pointers):
    //! - tree
    //! - sequence alignment
    //! - a BranchSelector of branch lengths
    //! - a site Selector of site rates (if pointer is null, then rates across
    //! sites are all equal to 1)
    //! - a PolyProcess to compute the likelihood of the data taking into account
    //! occurrences in the population of the reference and derived alleles
    PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
        const BranchSelector<double> *inbranchlength, const Selector<double> *insiterate,
        PolyProcess *inpolyprocess);


    //! \brief generic constructor
    //!
    //! Constructor takes as parameters (pointers):
    //! - tree
    //! - sequence alignment
    //! - a BranchSelector of branch lengths
    //! - a site Selector of site rates (if pointer is null, then rates across
    //! sites are all equal to 1)
    //! - a BranchSiteSelector specifying which substitution matrix should be used
    //! for each branch site pair
    //! - a site Selector of substitution matrices, specifying which matrix should
    //! be used for getting the equilibrium frequencies at each site, from which
    //! to draw the root state
    //! - a PolyProcess to compute the likelihood of the data taking into account
    //! occurrences in the population of the reference and derived alleles,
    //! this is nullpointer if no data could be found
    PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
        const BranchSelector<double> *inbranchlength, const Selector<double> *insiterate,
        const BranchSiteSelector<SubMatrix> *insubmatrixarray,
        const Selector<SubMatrix> *inrootsubmatrixarray, PolyProcess *inpolyprocess = nullptr);

    //! \brief special (short-cut) constructor for branch-homogeneous and
    //! site-homogeneous model
    //!
    //! Compared to the generic constructor, this constructor takes a pointer to a
    //! single substitution matrix. The branch/site and site selectors of
    //! substitution matrices (submatrixarray and rootsubmatrixarray) are then
    //! internally allocated by PhyloProcess based on this matrix. If insiterate
    //! pointer is null, then rates across sites are all equal to 1.
    //! - a PolyProcess to compute the likelihood of the data taking into account
    //! occurrences in the population of the reference and derived alleles,
    //! this is nullpointer if no data could be found
    PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
        const BranchSelector<double> *inbranchlength, const Selector<double> *insiterate,
        const SubMatrix *insubmatrix, PolyProcess *inpolyprocess = nullptr);

    //! \brief special (short-cut) constructor for branch-homogeneous and
    //! site-heterogeneous model
    //!
    //! Compared to the generic constructor, this constructor takes a pointer to a
    //! (site) Selector<SubMatrix>* insubmatrixarray. The branch/site selector of
    //! substitution matrices (submatrixarray) is then internally allocated by
    //! PhyloProcess based on this matrix based on this array, while
    //! rootmatrixarray is set to insubmatrixarray. If insiterate pointer is null,
    //! then rates across sites are all equal to 1.
    //! - a PolyProcess to compute the likelihood of the data taking into account
    //! occurrences in the population of the reference and derived alleles,
    //! this is nullpointer if no data could be found
    PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
        const BranchSelector<double> *inbranchlength, const Selector<double> *insiterate,
        const Selector<SubMatrix> *insubmatrixarray, PolyProcess *inpolyprocess = nullptr);

    //! \brief special (short-cut) constructor for branch-heterogeneous and
    //! site-homogeneous model
    //!
    //! Compared to the generic constructor, this constructor takes a pointer to
    //! BranchSelector<SubMatrix>* insubmatrixarray and a single submatrix (for
    //! the root eq frequencies). The branch/site and site selectors of
    //! substitution matrices (submatrixarray and rootsubmatrixarray) are then
    //! internally allocated by PhyloProcess based on these parameters. If
    //! insiterate pointer is null, then rates across sites are all equal to 1.
    //! - a PolyProcess to compute the likelihood of the data taking into account
    //! occurrences in the population of the reference and derived alleles,
    //! this is nullpointer if no data could be found
    PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
        const BranchSelector<double> *inbranchlength, const Selector<double> *insiterate,
        const BranchSelector<SubMatrix> *insubmatrixbrancharray, const SubMatrix *inrootsubmatrix,
        PolyProcess *inpolyprocess = nullptr);

    ~PhyloProcess();

    //! return log likelihood (computed using the pruning algorithm, Felsenstein
    //! 1981)
    double GetLogLikelihood() const;

    //! return log likelihood for given site
    double SiteLogLikelihood(int site) const;

    //! stochastic sampling of substitution history under current parameter
    //! configuration
    void ResampleSub();

    //! stochastic sampling of substitution history for given random fraction of
    //! sites
    double Move(double fraction);

    //! create all data structures necessary for computation
    void Unfold();

    //! delete data structures
    void Cleanup();

    //! posterior predictive resampling under current parameter configuration
    void PostPredSample(std::string name, bool rootprior = true);  // unclamped Nielsen

    //! get data from tips (after simulation) and put in into sequence alignment
    void GetLeafData(SequenceAlignment *data);

    int GetPathState(int taxon, int site) const {
        int node = reverse_taxon_table[taxon];
        auto site_leaf_path_map = pathmap[node][site];
        return site_leaf_path_map->GetFinalState();
    }

  private:
    double GetFastLogProb() const;
    double FastSiteLogLikelihood(int site) const;

    //! return branch length for given branch, based on index of node at the tip of the branch
    double GetBranchLength(Tree::NodeIndex index) const {
        return branchlength->GetVal(tree->branch_index(index));
    }

    //! return site rate for given site (if no rates-across-sites array was given
    //! to phyloprocess, returns 1)
    double GetSiteRate(int site) const {
        if (!siterate) { return 1.0; }
        return siterate->GetVal(site);
    }

    //! return matrix that should be used on a given branch based on index of node at branch tip
    const SubMatrix &GetSubMatrix(Tree::NodeIndex index, int site) const {
        return submatrixarray->GetVal(tree->branch_index(index), site);
    }

    const EVector &GetRootFreq(int site) const {
        return rootsubmatrixarray->GetVal(site).GetStationary();
    }

    const StateSpace *GetStateSpace() const { return data->GetStateSpace(); }
    const TaxonSet *GetTaxonSet() const { return data->GetTaxonSet(); }

    int GetNsite() const { return data->GetNsite(); }
    int GetNtaxa() const { return data->GetNtaxa(); }
    int GetNnode() const { return tree->nb_nodes(); }

    int GetNstate() const { return Nstate; }

    const SequenceAlignment *GetData() const { return data; }
    int GetData(int taxon, int site) const {
        if (taxon_table[taxon] == -1) {
            std::cerr << "error in taxon correspondance table\n";
            exit(1);
        }
        return data->GetState(taxon_table[taxon], site);
    }

    const Tree *GetTree() const { return tree; }
    Tree::NodeIndex GetRoot() const { return GetTree()->root(); }

    int GetMaxTrial() const { return maxtrial; }
    void SetMaxTrial(int i) { maxtrial = i; }

    void SetData(const SequenceAlignment *indata);
    void ClampData() { clampdata = true; }
    void UnclampData() { clampdata = false; }

    bool isDataCompatible(int taxon, int site, int state) const {
        return GetStateSpace()->isCompatible(GetData(taxon, site), state);
    }

    void DrawSites(double fraction);  // draw a fraction of sites which will be resampled
    void ResampleSub(int site);

    //! compute path sufficient statistics across all sites and branches and add
    //! them to suffstat (site-branch-homogeneous model)
    void AddPathSuffStat(PathSuffStat &suffstat) const;
    //! compute path sufficient statistics across all sites and branches and add
    //! them to suffstatarray (site-heterogeneous branch-homogeneous model)
    void AddPathSuffStat(Array<PathSuffStat> &suffstatarray) const;
    //! compute path sufficient statistics across all sites and branches and add
    //! them to suffstatarray (branch-heterogeneous site-homogeneous model)
    void AddPathSuffStat(NodeArray<PathSuffStat> &suffstatarray) const;
    //! compute path sufficient statistics across all sites and branches and add
    //! them to bidim suffstatarray (branches partitioned into conditions, see
    //! DiffSelModel) heterogeneeous across sites, branches partitioned into
    //! conditions
    void AddPathSuffStat(
        BidimArray<PathSuffStat> &suffstatarray, const BranchSelector<int> &branchalloc) const;

    //! compute path sufficient statistics across all sites and branches and add
    //! them to suffstat (site-branch-homogeneous model)
    void AddPolySuffStat(PolySuffStat &suffstat) const;
    //! compute path sufficient statistics across all sites and branches and add
    //! them to suffstatarray (site-heterogeneous branch-homogeneous model)
    void AddPolySuffStat(Array<PolySuffStat> &suffstatarray) const;

    //! compute path sufficient statistics for resampling branch lengths add them
    //! to branchlengthpathsuffstatarray
    void AddLengthSuffStat(BranchArray<PoissonSuffStat> &branchlengthpathsuffstatarray) const;

    //! compute path sufficient statistics for resampling branch lengths add them
    //! to branchlengthpathsuffstatarray
    void AddRateSuffStat(Array<PoissonSuffStat> &siteratepathsuffstatarray) const;

    void RecursiveAddPathSuffStat(Tree::NodeIndex from, PathSuffStat &suffstat) const;
    void LocalAddPathSuffStat(Tree::NodeIndex from, PathSuffStat &suffstat) const;

    void RecursiveAddPathSuffStat(
        Tree::NodeIndex from, NodeArray<PathSuffStat> &suffstatarray) const;
    void LocalAddPathSuffStat(Tree::NodeIndex from, NodeArray<PathSuffStat> &suffstatarray) const;

    void RecursiveAddPathSuffStat(Tree::NodeIndex from, Array<PathSuffStat> &suffstatarray) const;
    void LocalAddPathSuffStat(Tree::NodeIndex from, Array<PathSuffStat> &suffstatarray) const;

    void RecursiveAddPathSuffStat(Tree::NodeIndex from, BidimArray<PathSuffStat> &suffstatarray,
        const BranchSelector<int> &branchalloc) const;
    void LocalAddPathSuffStat(
        Tree::NodeIndex from, BidimArray<PathSuffStat> &suffstatarray, int cond) const;

    void RecursiveAddLengthSuffStat(
        Tree::NodeIndex from, BranchArray<PoissonSuffStat> &branchlengthpathsuffstatarray) const;
    void LocalAddLengthSuffStat(Tree::NodeIndex from, PoissonSuffStat &branchlengthsuffstat) const;

    void RecursiveAddRateSuffStat(
        Tree::NodeIndex from, Array<PoissonSuffStat> &siteratepathsuffstatarray) const;
    void LocalAddRateSuffStat(
        Tree::NodeIndex from, Array<PoissonSuffStat> &siteratepathsuffstatarray) const;

    void PostPredSample(int site, bool rootprior = false);
    // rootprior == true : root state drawn from stationary probability of the
    // process
    // rootprior == false: root state drawn from posterior distribution

    // various accessors

    /*
    bool isMissing(Tree::NodeIndex node, int site) const { return false; }
    bool isMissing(const Link *link, int site) const {
        return false;
        // return (missingmap[link->GetNode()][site] ||
        // missingmap[link->Out()->GetNode()][site]);
    }
    */

    void CreateMissingMap();
    void DeleteMissingMap();
    void RecursiveCreateMissingMap(Tree::NodeIndex from);
    void FillMissingMap();
    void BackwardFillMissingMap(Tree::NodeIndex from);
    void ForwardFillMissingMap(Tree::NodeIndex from, Tree::NodeIndex up);

    double GetPruningTime() const { return pruningchrono.GetTime(); }
    double GetResampleTime() const { return resamplechrono.GetTime(); }

    void RecursiveCreate(Tree::NodeIndex from);
    void RecursiveDelete(Tree::NodeIndex from);

    void RecursiveCreateTBL(Tree::NodeIndex from);
    void RecursiveDeleteTBL(Tree::NodeIndex from);

    void Pruning(Tree::NodeIndex from, int site) const;
    void ResampleSub(Tree::NodeIndex from, int site);
    void ResampleState();
    void ResampleState(int site);
    void PruningAncestral(Tree::NodeIndex from, int site);
    void PriorSample(Tree::NodeIndex from, int site, bool rootprior);
    void PriorSample();
    void RootPosteriorDraw(int site);

    // borrowed from phylobayes
    // where should that be?
    BranchSitePath *SamplePath(
        int stateup, int statedown, double time, double rate, const SubMatrix &matrix);
    BranchSitePath *SampleRootPath(int rootstate);
    BranchSitePath *ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate,
        double totaltime, const SubMatrix &matrix);
    BranchSitePath *ResampleUniformized(
        int stateup, int statedown, double rate, double totaltime, const SubMatrix &matrix);

    const Tree *tree;
    const SequenceAlignment *data;
    std::vector<int> taxon_table;
    std::vector<int> reverse_taxon_table;
    PolyProcess *polyprocess;
    const BranchSelector<double> *branchlength;
    const Selector<double> *siterate;
    const BranchSiteSelector<SubMatrix> *submatrixarray;
    const Selector<SubMatrix> *rootsubmatrixarray;
    bool allocsubmatrixarray;
    bool allocrootsubmatrixarray;

    int *sitearray;
    mutable double *sitelnL;

    int Nstate;

    bool clampdata;

    mutable double **uppercondlmap;
    mutable double **lowercondlmap;
    mutable BranchSitePath ***pathmap;
    int **statemap;
    int **missingmap;

    int maxtrial;
    static const int unknown = -1;

    static const int DEFAULTMAXTRIAL = 100;

    mutable Chrono pruningchrono;
    mutable Chrono resamplechrono;
};