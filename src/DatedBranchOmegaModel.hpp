#pragma once

#include "BranchProduct.hpp"
#include "Chronogram.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixBranchArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "Move.hpp"
#include "NodeProcess.hpp"
#include "PhyloProcess.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"

/**
 * \brief A site-homogeneous and branch-heterogenous Muse and Gaut omega-codon model.
 * Mutation rate and omega are modeled as a browian process along the tree (like CoEvol)
 *
 * The model has the following structure:
 * - nucleotide relative exchangeabilities and stationaries are uniform
 * Dirichlet
 * - there is one omega=dN/dS for each branches but for all sites
 * - branch rates (mutation rates) and omega are brownian multivariate
 */


class DatedBranchOmegaModel : public ChainComponent {
    // tree and data
    std::string datafile, treefile;
    std::unique_ptr<const Tree> tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    // Node ages
    NodeAges *nodeages;
    // Chronogram (diff between node ages)
    Chronogram *chronogram;

    // Branch rates (brownian process)
    double ratesprocesssigma2;
    double ratesprocessroot;
    NodeProcess *nodeprocessrates;
    BranchProcess *branchprocessrates;

    // Branch lengths (product of branch rates and chronogram)
    BranchwiseProduct<double> *branchlength;

    // Branch omega (brownian process)
    double omegaprocesssigma2;
    double omegaprocessroot;
    NodeProcess *nodeprocessomega;
    BranchProcess *branchprocessomega;

    // Nucleotide rates
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucrelrate;
    std::vector<double> nucstat;
    GTRSubMatrix *nucmatrix;

    // codon matrices (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixBranchArray *codonmatrixbrancharray;
    MGOmegaCodonSubMatrix *rootcodonmatrix;

    // PhyloProcess
    PhyloProcess *phyloprocess;

    // suff stats for substitution paths
    // summed over all branches and over all sites
    PathSuffStatNodeArray *pathsuffstatarray;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function
    // of nodeprocessomega
    OmegaPathSuffStatBranchArray *omegapathsuffstatarray;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

  public:
    friend std::ostream &operator<<(std::ostream &os, DatedBranchOmegaModel &m);

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    DatedBranchOmegaModel(std::string datafile, std::string treefile)
        : datafile(datafile), treefile(treefile) {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);
        Nbranch = tree->nb_branches();

        // Node ages
        nodeages = new NodeAges(*tree);
        // Chronogram (diff between node ages)
        chronogram = new Chronogram(*nodeages);

        // Branch rates (brownian process)
        ratesprocesssigma2 = 1.0;
        ratesprocessroot = 0.1;
        nodeprocessrates = new NodeProcess(*chronogram, ratesprocesssigma2, ratesprocessroot);
        branchprocessrates = new BranchProcess(*nodeprocessrates);

        // Branch lengths (product of branch rates and chronogram)
        branchlength = new BranchwiseProduct<double>(*chronogram, *branchprocessrates);

        // Branch omega (brownian process)
        omegaprocessroot = 0.3;
        omegaprocesssigma2 = 1.0;
        nodeprocessomega = new NodeProcess(*chronogram, omegaprocesssigma2, omegaprocessroot);
        branchprocessomega = new BranchProcess(*nodeprocessomega);

        // Nucleotide rates
        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        // Codon Matrices
        codonmatrixbrancharray = new MGOmegaCodonSubMatrixBranchArray(
            GetCodonStateSpace(), nucmatrix, branchprocessomega);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(
            GetCodonStateSpace(), nucmatrix, nodeprocessomega->GetVal(tree->root()));

        // PhyloProcess
        phyloprocess = new PhyloProcess(
            tree.get(), codondata, branchlength, nullptr, codonmatrixbrancharray, rootcodonmatrix);
        phyloprocess->Unfold();

        // Suff Stats
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        omegapathsuffstatarray = new OmegaPathSuffStatBranchArray(*tree);
        pathsuffstatarray = new PathSuffStatNodeArray(*tree);
    }

    virtual ~DatedBranchOmegaModel() = default;

    void move(int it) override { Move(); }

    template <class C>
    void declare_model(C &t) {
        t.add("nucstat", nucstat);
        t.add("nucrelrate", nucrelrate);
        //t.add("nodeages", *nodeages);
        //t.add("nodeprocessrates", *nodeprocessrates);
        //t.add("nodeprocessomega", *nodeprocessomega);
        t.add("ratesprocessroot", ratesprocessroot);
        t.add("ratesprocesssigma2", ratesprocesssigma2);
        t.add("omegaprocessroot", omegaprocessroot);
        t.add("omegaprocesssigma2", omegaprocesssigma2);
    }

    template <class C>
    void declare_stats(C &t) {
        t.add("logprior", this, &DatedBranchOmegaModel::GetLogPrior);
        t.add("lnL", this, &DatedBranchOmegaModel::GetLogLikelihood);
        t.add("bomegamean", [&]() { return branchprocessomega->GetMean(); });
        t.add("bratesmean", [&]() { return branchprocessomega->GetVar(); });
        t.add("bratesmean", [&]() { return branchprocessrates->GetMean(); });
        t.add("bratesvar", [&]() { return branchprocessrates->GetVar(); });
        t.add("statent", [&]() { return Random::GetEntropy(nucstat); });
        t.add("rrent", [&]() { return Random::GetEntropy(nucrelrate); });
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrix that its parameters have changed and that it
    //! should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchCodonMatrix() {
        codonmatrixbrancharray->UpdateCodonMatrices();
        rootcodonmatrix->SetOmega(nodeprocessomega->GetVal(tree->root()));
        rootcodonmatrix->CorruptMatrix();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices() {
        TouchNucMatrix();
        TouchCodonMatrix();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in Move), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief Update the BranchProcess (rates and omega) with the underlying NodeProcess, Update
    //! the Chronogram with the underlying NodeAges. And finally update the branch lengths with the
    //! Chronogram and the BranchProcess (rates).
    //!
    //! Used when the model is restarted or for the posterior predictif.
    void UpdateBranches() {
        chronogram->Update();
        branchprocessomega->Update();
        branchprocessrates->Update();
        branchlength->Update();
    }

    //! \brief Update the chronogram (branch time) and branch lengths around the focal node.
    //!
    //! Update needed when the age (NodeAges) of the focal node is changed.
    void UpdateLocalChronogram(Tree::NodeIndex node) {
        chronogram->UpdateLocal(node);
        branchlength->UpdateLocal(node);
    }

    //! \brief Update the branch rates and lengths around the focal node.
    //!
    //! Update needed when the rate (NodeProcess) of the focal node is changed.
    void UpdateLocalBranchProcessRates(Tree::NodeIndex node) {
        branchprocessrates->UpdateLocal(node);
        branchlength->UpdateLocal(node);
    }

    //! \brief Update the branch omega around the focal node.
    //!
    //! Update needed when the omega (NodeProcess) of the focal node is changed.
    void UpdateLocalBranchProcessOmega(Tree::NodeIndex node) {
        branchprocessomega->UpdateLocal(node);
        if (!tree->is_root(node)) {
            Tree::BranchIndex branch = tree->branch_index(node);
            (*codonmatrixbrancharray)[branch].SetOmega(branchprocessomega->GetVal(branch));
            (*codonmatrixbrancharray)[branch].CorruptMatrix();
        }
    }

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() {
        UpdateBranches();
        TouchMatrices();
        ResampleSub(1.0);
    }

    //-------------------
    // Posterior Predictive
    // ------------------

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(std::string name) {
        UpdateBranches();
        TouchMatrices();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        total += NodeProcessRatesLogPrior();
        total += NodeProcessRatesHyperLogPrior();
        total += NodeProcessOmegaLogPrior();
        total += NodeProcessOmegaHyperLogPrior();
        total += NucRatesLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    // branche rates prior

    //! log prior over branch rates (brownian process)
    double NodeProcessRatesLogPrior() const { return nodeprocessrates->GetLogProb(); }

    //! log prior of branch rate (brownian process) around of focal node
    double LocalNodeProcessRatesLogPrior(Tree::NodeIndex node) const {
        return nodeprocessrates->GetLocalLogProb(node);
    }

    //! log prior of hyperparameters of branch rates (brownian process)
    double NodeProcessRatesHyperLogPrior() const {
        double total = 0;
        // exponential of mean 1 on shape
        total -= ratesprocessroot;
        // exponential of mean 1 on rate
        total += Random::logInverseGammaDensity(ratesprocesssigma2, 2, 1.0 / 2);
        return total;
    }

    //! log prior over branch omega (brownian process)
    double NodeProcessOmegaLogPrior() const { return nodeprocessomega->GetLogProb(); }

    //! log prior of branch omega (brownian process) around of focal node
    double LocalNodeProcessOmegaLogPrior(Tree::NodeIndex node) const {
        return nodeprocessomega->GetLocalLogProb(node);
    }

    //! log prior of hyperparameters of branch omega (brownian process)
    double NodeProcessOmegaHyperLogPrior() const {
        double total = 0;
        // exponential of mean 1 on shape
        total -= omegaprocessroot;
        // exponential of mean 1 on rate
        total += Random::logInverseGammaDensity(omegaprocesssigma2, 2, 1.0 / 2);
        return total;
    }

    // Nucleotide rates prior
    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq.
    //! freqs. (nucstat) -- uniform Dirichlet in both cases
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(
            nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //-------------------
    // Collect Suff Stat
    //-------------------

    // Paths
    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    // Branch lengths
    //! collect sufficient statistics for moving branch lengths
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    // Nucleotide rates
    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixbrancharray, *rootcodonmatrix, *pathsuffstatarray);
    }

    // Omega
    //! collect sufficient statistics for moving omega (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectOmegaPathSuffStat() {
        omegapathsuffstatarray->Clear();
        omegapathsuffstatarray->AddSuffStat(
            *codonmatrixbrancharray, *rootcodonmatrix, *pathsuffstatarray);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // Node ages and branch rates
    //! \brief log prob to be recomputed when moving age of focal node
    double LocalNodeAgeLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        tot += LocalNodeProcessRatesLogPrior(node);
        tot += LocalNodeProcessOmegaLogPrior(node);
        tot += LocalBranchLengthSuffStatLogProb(node);
        return tot;
    }

    //! \brief log prob to be recomputed when moving branch rates (brownian process) around of focal
    //! node
    double LocalNodeProcessRatesLogProb(Tree::NodeIndex node) const {
        return LocalNodeProcessRatesLogPrior(node) + LocalBranchLengthSuffStatLogProb(node);
    }

    //! \brief log prob factor (without prior) to be recomputed when moving age of focal node, or
    //! when moving branch rates (brownian process) around of focal node.
    double LocalBranchLengthSuffStatLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        // for all children
        for (auto const &child : tree->children(node)) {
            tot += BranchLengthSuffStatLogProb(tree->branch_index(child));
        }
        if (!tree->is_root(node)) {
            // for the branch attached to the node
            tot += BranchLengthSuffStatLogProb(tree->branch_index(node));
        }
        assert(tot != 0);
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! the length of a given branch
    double BranchLengthSuffStatLogProb(Tree::BranchIndex branch) const {
        return lengthpathsuffstatarray->GetVal(branch).GetLogProb(branchlength->GetVal(branch));
    }

    // Omega
    //! \brief log prob to be recomputed when moving omega (brownian process) around of focal node
    double LocalNodeProcessOmegaLogProb(Tree::NodeIndex node) const {
        return LocalNodeProcessOmegaLogPrior(node) + LocalNodeProcessOmegaSuffStatLogProb(node);
    }

    //! \brief log prob factor (without prior) to be recomputed when moving omega (brownian process)
    //! around of focal node
    double LocalNodeProcessOmegaSuffStatLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        // for all children
        for (auto const &child : tree->children(node)) {
            tot += NodeProcessOmegaSuffStatLogProb(tree->branch_index(child));
        }
        if (!tree->is_root(node)) {
            // for the branch attached to the node
            tot += NodeProcessOmegaSuffStatLogProb(tree->branch_index(node));
        }
        assert(tot != 0);
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! omega of a given branch
    double NodeProcessOmegaSuffStatLogProb(Tree::BranchIndex branch) const {
        return omegapathsuffstatarray->GetVal(branch).GetLogProb(
            branchprocessomega->GetVal(branch));
    }


    //! \brief log prob factor to be recomputed when moving branch rates (brownian process)
    //! hyperparams
    double NodeProcessRatesHyperLogProb() const {
        return NodeProcessRatesHyperLogPrior() + NodeProcessRatesLogPrior();
    }

    //! \brief log prob factor to be recomputed when moving omega (brownian process) hyperparams
    double NodeProcessOmegaHyperLogProb() const {
        return NodeProcessOmegaHyperLogPrior() + NodeProcessOmegaLogPrior();
    }

    // Nucleotide rates
    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    //! Calculated using nucpathsuffstat
    //! (which summarizes all information about how the probability of the
    //! substitution mapping depends on nucleotide mutation rates) and the
    //! nucmatrix. Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        TouchMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            CollectLengthSuffStat();
            MoveNodeAges(1.0, 3);
            MoveNodeAges(0.1, 3);

            MoveNodeProcessRates(1.0, 3);
            MoveNodeProcessRates(0.1, 3);
            MoveNodeProcessRatesHyperParameters();

            CollectPathSuffStat();
            CollectNucPathSuffStat();
            MoveNucRates();
            TouchMatrices();

            CollectOmegaPathSuffStat();
            MoveNodeProcessOmega(1.0, 3);
            MoveNodeProcessOmega(0.1, 3);
            MoveNodeProcessOmegaHyperParameters();
            TouchMatrices();
        }
    }

    //! MH moves on branch ages
    void MoveNodeAges(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(tree->nb_nodes()); node++) {
                if (!tree->is_root(node) and !tree->is_leaf(node)) { MoveNodeAge(node, tuning); }
            }
        }
    }

    //! MH moves on branch ages for a focal node
    void MoveNodeAge(Tree::NodeIndex node, double tuning) {
        double bk = nodeages->GetVal(node);
        double logratio = -LocalNodeAgeLogProb(node);
        double sliding = tuning * (Random::Uniform() - 0.5);

        nodeages->SlidingMove(node, sliding);
        UpdateLocalChronogram(node);

        logratio += LocalNodeAgeLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            (*nodeages)[node] = bk;
            UpdateLocalChronogram(node);
        }
    }

    //! MH moves on branch rates (brownian process)
    void MoveNodeProcessRates(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(tree->nb_nodes()); node++) {
                MoveNodeProcessRate(node, tuning);
            }
        }
    }

    //! MH moves on branch rates (brownian process) for a focal node
    void MoveNodeProcessRate(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodeProcessRatesLogProb(node);

        double m = tuning * sqrt(ratesprocesssigma2) * (Random::Uniform() - 0.5);
        (*nodeprocessrates)[node] += m;
        UpdateLocalBranchProcessRates(node);

        logratio += LocalNodeProcessRatesLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            (*nodeprocessrates)[node] -= m;
            UpdateLocalBranchProcessRates(node);
        }
    }

    //! MH moves on branch rates (brownian process) hyperparameters
    void MoveNodeProcessRatesHyperParameters() {
        Move::Scaling(ratesprocesssigma2, 1.0, 3,
            &DatedBranchOmegaModel::NodeProcessRatesHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
        Move::Scaling(ratesprocesssigma2, 0.1, 3,
            &DatedBranchOmegaModel::NodeProcessRatesHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
        Move::Scaling(ratesprocessroot, 1.0, 3,
            &DatedBranchOmegaModel::NodeProcessRatesHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
        Move::Scaling(ratesprocessroot, 0.1, 3,
            &DatedBranchOmegaModel::NodeProcessRatesHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
    }

    //! MH moves on branch omega (brownian process)
    void MoveNodeProcessOmega(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(tree->nb_nodes()); node++) {
                MoveNodeOmega(node, tuning);
            }
        }
    }

    //! MH moves on branch omega (brownian process) for a focal node
    void MoveNodeOmega(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodeProcessOmegaLogProb(node);

        double m = tuning * sqrt(omegaprocesssigma2) * (Random::Uniform() - 0.5);
        (*nodeprocessomega)[node] += m;

        UpdateLocalBranchProcessOmega(node);
        logratio += LocalNodeProcessOmegaLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            (*nodeprocessomega)[node] -= m;
            UpdateLocalBranchProcessOmega(node);
        }
    }

    //! MH moves on branch omega (brownian process) hyperparameters
    void MoveNodeProcessOmegaHyperParameters() {
        Move::Scaling(omegaprocesssigma2, 1.0, 3,
            &DatedBranchOmegaModel::NodeProcessOmegaHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
        Move::Scaling(omegaprocesssigma2, 1.0, 3,
            &DatedBranchOmegaModel::NodeProcessOmegaHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
        Move::Scaling(omegaprocessroot, 1.0, 3,
            &DatedBranchOmegaModel::NodeProcessOmegaHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
        Move::Scaling(omegaprocessroot, 0.1, 3,
            &DatedBranchOmegaModel::NodeProcessOmegaHyperLogProb, &DatedBranchOmegaModel::NoUpdate,
            this);
    }

    // Nucleotide rates
    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        Move::Profile(nucrelrate, 0.1, 1, 3, &DatedBranchOmegaModel::NucRatesLogProb,
            &DatedBranchOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucrelrate, 0.03, 3, 3, &DatedBranchOmegaModel::NucRatesLogProb,
            &DatedBranchOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucrelrate, 0.01, 3, 3, &DatedBranchOmegaModel::NucRatesLogProb,
            &DatedBranchOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucstat, 0.1, 1, 3, &DatedBranchOmegaModel::NucRatesLogProb,
            &DatedBranchOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucstat, 0.01, 1, 3, &DatedBranchOmegaModel::NucRatesLogProb,
            &DatedBranchOmegaModel::TouchNucMatrix, this);
    }

    void ToStream(std::ostream &os) { os << *this; }
};

std::istream &operator>>(std::istream &is, std::unique_ptr<DatedBranchOmegaModel> &m) {
    std::string model_name;
    std::string datafile;
    std::string treefile;

    is >> model_name;
    if (model_name != "DatedBranchOmega") {
        std::cerr << "Expected DatedBranchOmega for model name, got " << model_name << "\n";
        exit(1);
    }

    is >> datafile;
    is >> treefile;
    m.reset(new DatedBranchOmegaModel(datafile, treefile));
    Tracer tracer{*m, &DatedBranchOmegaModel::declare_model};
    tracer.read_line(is);
    m->Update();
    return is;
}

std::ostream &operator<<(std::ostream &os, DatedBranchOmegaModel &m) {
    Tracer tracer{m, &DatedBranchOmegaModel::declare_model};
    os << "DatedBranchOmega" << '\t';
    os << m.datafile << '\t';
    os << m.treefile << '\t';
    tracer.write_line(os);
    return os;
}
