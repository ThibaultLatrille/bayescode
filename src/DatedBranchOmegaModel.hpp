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
#include "NodeMultivariateProcess.hpp"
#include "PhyloProcess.hpp"
#include "ScatterSuffStat.hpp"
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

    int dimension;
    int cov_df;
    double cov_kappa;
    // Covariance matrix
    CovMatrix *cov_matrix;
    EVector *root_multivariate;
    NodeMultivariateProcess *node_multivariate;

    // Branch rates (brownian process)
    NodeProcess *noderates;
    BranchProcess *branchrates;

    // Branch lengths (product of branch rates and chronogram)
    BranchwiseProduct<double> *branchlength;

    // Branch omega (brownian process)
    NodeProcess *nodeomega;
    BranchProcess *branchomega;

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
    // of nodeomega
    OmegaPathSuffStatBranchArray *omegapathsuffstatarray;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    ScatterSuffStat *scattersuffstat;

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

        dimension = 2;
        cov_df = dimension + 1;
        cov_kappa = 1.0;
        root_multivariate = new EVector(dimension);
        *root_multivariate = EVector::Constant(dimension, 0.1);
        cov_matrix = new CovMatrix(dimension);

        node_multivariate = new NodeMultivariateProcess(*chronogram, *cov_matrix, *root_multivariate);

        // Branch rates (brownian process)
        noderates = new NodeProcess(*node_multivariate, 0);
        branchrates = new BranchProcess(*noderates);

        // Branch lengths (product of branch rates and chronogram)
        branchlength = new BranchwiseProduct<double>(*chronogram, *branchrates);

        // Branch omega (brownian process)
        nodeomega = new NodeProcess(*node_multivariate, 1);
        branchomega = new BranchProcess(*nodeomega);

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
        codonmatrixbrancharray =
            new MGOmegaCodonSubMatrixBranchArray(GetCodonStateSpace(), nucmatrix, branchomega);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(
            GetCodonStateSpace(), nucmatrix, nodeomega->GetVal(tree->root()));

        // PhyloProcess
        phyloprocess = new PhyloProcess(
            tree.get(), codondata, branchlength, nullptr, codonmatrixbrancharray, rootcodonmatrix);
        phyloprocess->Unfold();

        // Suff Stats
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        omegapathsuffstatarray = new OmegaPathSuffStatBranchArray(*tree);
        pathsuffstatarray = new PathSuffStatNodeArray(*tree);
        scattersuffstat = new ScatterSuffStat(*tree);
    }

    virtual ~DatedBranchOmegaModel() = default;

    void move(int it) override { Move(); }

    template <class C>
    void declare_model(C &t) {
        t.add("nucstat", nucstat);
        t.add("nucrelrate", nucrelrate);
        // t.add("nodeages", *nodeages);
        // t.add("covmatrix", *cov_matrix);
        // t.add("root_multivariate", root_multivariate);
        // t.add("node_multivariate", node_multivariate);
        t.add("cov_kappa", cov_kappa);
        t.add("cov_df", cov_df);
    }

    template <class C>
    void declare_stats(C &t) {
        t.add("logprior", this, &DatedBranchOmegaModel::GetLogPrior);
        t.add("lnL", this, &DatedBranchOmegaModel::GetLogLikelihood);
        t.add("bomegamean", [&]() { return branchomega->GetMean(); });
        t.add("bratesmean", [&]() { return branchomega->GetVar(); });
        t.add("bratesmean", [&]() { return branchrates->GetMean(); });
        t.add("bratesvar", [&]() { return branchrates->GetVar(); });
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
        rootcodonmatrix->SetOmega(nodeomega->GetVal(tree->root()));
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
        branchomega->Update();
        branchrates->Update();
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
    void UpdateLocalBranchRates(Tree::NodeIndex node) {
        branchrates->UpdateLocal(node);
        branchlength->UpdateLocal(node);
    }

    //! \brief Update the branch omega around the focal node.
    //!
    //! Update needed when the omega (NodeProcess) of the focal node is changed.
    void UpdateLocalBranchOmega(Tree::NodeIndex node) {
        branchomega->UpdateLocal(node);
        if (!tree->is_root(node)) {
            Tree::BranchIndex branch = tree->branch_index(node);
            (*codonmatrixbrancharray)[branch].SetOmega(branchomega->GetVal(branch));
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
        total += NodeMultivariateLogPrior();
        total += RootMultivariateLogPrior();
        total += NucRatesLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    // Multivariate prior

    //! log prior over branch rates (brownian process)
    double NodeMultivariateLogPrior() const { return node_multivariate->GetLogProb(); }

    //! log prior of branch rate (brownian process) around of focal node
    double LocalNodeMultivariateLogPrior(Tree::NodeIndex node) const {
        return node_multivariate->GetLocalLogProb(node);
    }

    //! log prior of
    double RootMultivariateLogPrior() const { return node_multivariate->GetLogProb(tree->root()); }

    //! log prior of hyperparameters of branch rates (brownian process)
    double RootMultivariateHyperLogPrior() const { return -root_multivariate->sum(); }

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

    // Scatter (brownian process)

    void CollectScatterSuffStat() {
        scattersuffstat->Clear();
        scattersuffstat->AddSuffStat(*node_multivariate);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // Node ages and branch rates
    //! \brief log prob to be recomputed when moving age of focal node
    double LocalNodeAgeLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        tot += LocalNodeMultivariateLogPrior(node);
        tot += LocalBranchLengthSuffStatLogProb(node);
        return tot;
    }

    //! \brief log prob to be recomputed when moving branch rates (brownian process) around of focal
    //! node
    double LocalNodeRatesLogProb(Tree::NodeIndex node) const {
        return LocalNodeMultivariateLogPrior(node) + LocalBranchLengthSuffStatLogProb(node);
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
    double LocalNodeOmegaLogProb(Tree::NodeIndex node) const {
        return LocalNodeMultivariateLogPrior(node) + LocalNodeOmegaSuffStatLogProb(node);
    }

    //! \brief log prob factor (without prior) to be recomputed when moving omega (brownian process)
    //! around of focal node
    double LocalNodeOmegaSuffStatLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        // for all children
        for (auto const &child : tree->children(node)) {
            tot += NodeOmegaSuffStatLogProb(tree->branch_index(child));
        }
        if (!tree->is_root(node)) {
            // for the branch attached to the node
            tot += NodeOmegaSuffStatLogProb(tree->branch_index(node));
        }
        assert(tot != 0);
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! omega of a given branch
    double NodeOmegaSuffStatLogProb(Tree::BranchIndex branch) const {
        return omegapathsuffstatarray->GetVal(branch).GetLogProb(branchomega->GetVal(branch));
    }

    //! Covariance Matrix log prob (log prior of the cov matrix + log prob of the nodeprocesses
    //! given the cov matrix
    double CovMatrixLogProb() const {
        return scattersuffstat->GetLogPosterior(*cov_matrix, cov_df, cov_kappa);
    }

    //! Root log prob
    double RootLogProb() const {
        return RootMultivariateHyperLogPrior() + RootMultivariateLogPrior();
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

            MoveNodeRates(1.0, 3);
            MoveNodeRates(0.1, 3);

            CollectPathSuffStat();
            CollectNucPathSuffStat();
            MoveNucRates();
            TouchMatrices();

            CollectOmegaPathSuffStat();
            MoveNodeOmega(1.0, 3);
            MoveNodeOmega(0.1, 3);

            CollectScatterSuffStat();
            //MoveCovMatrix(1.0, dimension * dimension);
            //MoveCovMatrix(0.1, dimension * dimension);
            MoveRootMultivariate(1.0, 3);
            MoveRootMultivariate(0.1, 3);
            TouchMatrices();
        }
    }

    void MoveCovMatrix(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            int i = Random::Choose(dimension);
            int j = Random::Choose(dimension);
            assert(i < dimension);
            assert(0 <= i);
            assert(j < dimension);
            assert(0 <= j);
            double bk = (*cov_matrix)(i, j);

            double logratio = -CovMatrixLogProb();
            double sliding = tuning * noderates->GetSigma() * (Random::Uniform() - 0.5);
            cov_matrix->SlidingMove(i, j, sliding);
            logratio += CovMatrixLogProb();

            bool accept = (log(Random::Uniform()) < logratio);
            if (!accept) { (*cov_matrix)(i, j) = bk; }
        }
    };

    void MoveRootMultivariate(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (int dim = 0; dim < dimension; dim++) {
                double logratio = -RootLogProb();

                double loghastings = tuning * (Random::Uniform() - 0.5);
                double hastings = exp(loghastings);
                (*root_multivariate)(dim) = hastings * (*root_multivariate)(dim);

                logratio += loghastings;
                logratio += RootLogProb();

                bool accept = (log(Random::Uniform()) < logratio);
                if (!accept) { (*root_multivariate)(dim) = (*root_multivariate)(dim) / hastings; }
            }
        }
    };

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
    void MoveNodeRates(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(tree->nb_nodes()); node++) {
                MoveNodeProcessRate(node, tuning);
            }
        }
    }

    //! MH moves on branch rates (brownian process) for a focal node
    void MoveNodeProcessRate(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodeRatesLogProb(node);

        double m = tuning * noderates->GetSigma() * (Random::Uniform() - 0.5);
        noderates->SlidingMove(node, m);
        UpdateLocalBranchRates(node);

        logratio += LocalNodeRatesLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            noderates->SlidingMove(node, m);
            UpdateLocalBranchRates(node);
        }
    }

    //! MH moves on branch omega (brownian process)
    void MoveNodeOmega(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(tree->nb_nodes()); node++) {
                MoveNodeOmega(node, tuning);
            }
        }
    }

    //! MH moves on branch omega (brownian process) for a focal node
    void MoveNodeOmega(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodeOmegaLogProb(node);

        double m = tuning * nodeomega->GetSigma() * (Random::Uniform() - 0.5);
        nodeomega->SlidingMove(node, m);

        UpdateLocalBranchOmega(node);
        logratio += LocalNodeOmegaLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            nodeomega->SlidingMove(node, -m);
            UpdateLocalBranchOmega(node);
        }
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
