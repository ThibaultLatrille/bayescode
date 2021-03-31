#pragma once

#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "Move.hpp"
#include "PhyloProcess.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"

/**
 * \brief A standard site- and branch-homogeneous Muse and Gaut omega-codon
 * model
 *
 * The model has the following structure:
 * - branch lengths iid Exponential of mean blhypermean
 * - nucleotide relative exchangeabilities and stationaries are uniform
 * Dirichlet
 * - there is one single omega=dN/dS for all sites and across all branches
 * - prior over omega ~ Gamma(omegahypermean,omegahyperinvshape).
 *
 * The 2 hyperparameters omegahypermean and hyperinvshape are fixed when this
 * model is used in isolation. On the other hand, this model can be used in a
 * multigene context (see MultiGeneSingeOmegaModel), in which case
 * omegahypermean and hyperinvshape are estimated across genes.
 */

// mode shared:      global
// mode shrunken:    gene specific, with hyperparameters estimated across genes
// mode independent: gene-specific, with fixed hyperparameters
enum param_mode_t { shared, shrunken, independent };

class SingleOmegaModel : public ChainComponent {
    // tree and data
    std::string datafile, treefile;
    std::unique_ptr<const Tree> tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    param_mode_t blmode;
    param_mode_t nucmode;

    // Branch lengths
    double blhypermean;
    double blhyperinvshape;
    SimpleBranchArray<double> *blhypermeanarray;
    GammaWhiteNoise *branchlength;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    // Nucleotide rates
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucrelrate;
    std::vector<double> nucstat;
    GTRSubMatrix *nucmatrix;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // Omega

    double omegahypermean;
    double omegahyperinvshape;
    double omega;

    // a codon matrix (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrix *codonmatrix;

    // PhyloProcess

    PhyloProcess *phyloprocess;

    // suff stats for substitution paths
    // summed over all branches and over all sites
    PathSuffStat pathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function
    // of omega
    OmegaPathSuffStat omegapathsuffstat;


  public:
    friend std::ostream &operator<<(std::ostream &os, SingleOmegaModel &m);

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    SingleOmegaModel(std::string datafile, std::string treefile, param_mode_t blmode = independent,
        param_mode_t nucmode = independent)
        : datafile(datafile), treefile(treefile), blmode(blmode), nucmode(nucmode) {
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

        // Branch lengths
        blhypermean = 0.1;
        blhyperinvshape = 1.0;
        blhypermeanarray = new SimpleBranchArray<double>(*tree, blhypermean);
        branchlength = new GammaWhiteNoise(*tree, *blhypermeanarray, blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);

        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        // Omega

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        omega = 1.0;
        codonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, omega);

        phyloprocess = new PhyloProcess(tree.get(), codondata, branchlength, 0, codonmatrix);
        phyloprocess->Unfold();
    }

    virtual ~SingleOmegaModel() = default;

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        model_node(info, "omega", omega);
        model_node(info, "nucstat", nucstat);
        model_node(info, "nucrelrate", nucrelrate);
        model_node(info, "branchlength", *branchlength);

        model_stat(info, "omega", omega);
        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        model_stat(info, "length", [this]() { return branchlength->GetTotalLength(); });
        model_stat(info, "statent", [&]() { return Random::GetEntropy(nucstat); });
        model_stat(info, "rrent", [&]() { return Random::GetEntropy(nucrelrate); });
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return current value of omega
    double GetOmega() const { return omega; }

    //-------------------
    // Setting and updating
    // ------------------

    // Branch lengths

    //! set branch lengths to a new value (multi-gene analyses)
    void SetBranchLengths(const BranchSelector<double> &inbranchlength) {
        branchlength->Copy(inbranchlength);
    }

    //! get a copy of branch lengths into array given as argument
    void GetBranchLengths(BranchArray<double> &inbranchlength) const {
        inbranchlength.Copy(*branchlength);
    }

    //! set branch lengths hyperparameters to a new value (multi-gene analyses)
    void SetBranchLengthsHyperParameters(
        const BranchSelector<double> &inblmeanarray, double inblinvshape) {
        blhypermeanarray->Copy(inblmeanarray);
        blhyperinvshape = inblinvshape;
    }

    // Nucleotide rates

    //! set nucleotide rates (relative exchangeabilities and eq. frequencies) to a
    //! new value (multi-gene analyses)
    void SetNucRates(
        const std::vector<double> &innucrelrate, const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        TouchMatrices();
    }

    //! get a copy of nucleotide rates into arrays given as arguments
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double> &innucrelratehypercenter,
        double innucrelratehyperinvconc, const std::vector<double> &innucstathypercenter,
        double innucstathyperinvconc) {
        nucrelratehypercenter = innucrelratehypercenter;
        nucrelratehyperinvconc = innucrelratehyperinvconc;
        nucstathypercenter = innucstathypercenter;
        nucstathyperinvconc = innucstathyperinvconc;
    }

    // Omega

    //! \brief set omega to a new value
    //!
    //! Used in a multigene context.
    //! Notifies corruption to the codon matrix.
    void SetOmega(double inomega) {
        omega = inomega;
        TouchCodonMatrix();
    }

    //! \brief set the hyperparameters of the gamma prior over omega
    //!
    //! Used in a multigene context.
    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape) {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
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
        codonmatrix->UpdateOmega(omega);
        codonmatrix->CorruptMatrix();
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
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() {
        TouchMatrices();
        ResampleSub(1.0);
    }

    //-------------------
    // Posterior Predictive
    // ------------------

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(std::string name) {
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
        if (blmode != shared) { total += BranchLengthsLogPrior(); }
        if (nucmode != shared) { total += NucRatesLogPrior(); }
        total += OmegaLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    // Branch lengths

    //! log prior over branch lengths
    double BranchLengthsLogPrior() const { return branchlength->GetLogProb(); }

    // Nucleotide rates

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

    // Omega

    //! log prior over omega (gamma of mean omegahypermean and shape
    //! 1/omegahyperinvshape)
    double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return Random::logGammaDensity(omega, alpha, beta);
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    // Branch lengths

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current
    //! substitution mapping
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    // Nucleotide rates

    //! \brief const acess to nuc-pathsuffstat
    //!
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate)
    //! and equilibrium frequencies (nucstat) conditional on the current
    //! substitution mapping.
    const NucPathSuffStat &GetNucPathSuffStat() const { return nucpathsuffstat; }

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrix, pathsuffstat);
    }

    //! \brief return log prob of current substitution mapping, as a function of
    //! nucleotide parameters (nucrelrate and nucstat)
    //!
    //! Calculated using nucpathsuffstat
    //! (which summarizes all information about how the probability of the
    //! substitution mapping depends on nucleotide mutation rates) and the
    //! nucmatrix. Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    // Paths

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstat.Clear();
        pathsuffstat.AddSuffStat(*phyloprocess);
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrix. Both pathsuffstat and
    //! codonmatrix are assumed to be updated.
    double PathSuffStatLogProb() const { return pathsuffstat.GetLogProb(*codonmatrix); }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // Nucleotide rates

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

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
            if (blmode != shared) { MoveBranchLengths(); }

            CollectPathSuffStat();
            MoveOmega();

            if (nucmode != shared) {
                TouchMatrices();
                MoveNucRates();
            }
        }
    }

    // Branch lengths

    //! overall schedule branch length updatdes
    void MoveBranchLengths() { ResampleBranchLengths(); }

    //! Gibbs resample branch lengths (based on sufficient statistics)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    // Nucleotide rates

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        Move::Profile(nucrelrate, 0.1, 1, 3, &SingleOmegaModel::NucRatesLogProb,
            &SingleOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucrelrate, 0.03, 3, 3, &SingleOmegaModel::NucRatesLogProb,
            &SingleOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucrelrate, 0.01, 3, 3, &SingleOmegaModel::NucRatesLogProb,
            &SingleOmegaModel::TouchNucMatrix, this);

        Move::Profile(nucstat, 0.1, 1, 3, &SingleOmegaModel::NucRatesLogProb,
            &SingleOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucstat, 0.01, 1, 3, &SingleOmegaModel::NucRatesLogProb,
            &SingleOmegaModel::TouchNucMatrix, this);

        TouchMatrices();
    }

    // Omega

    //! Gibbs resample omega (based on sufficient statistics of current
    //! substitution mapping)
    void MoveOmega() {
        omegapathsuffstat.Clear();
        omegapathsuffstat.AddSuffStat(*codonmatrix, pathsuffstat);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(
            alpha + omegapathsuffstat.GetCount(), beta + omegapathsuffstat.GetBeta());
        TouchCodonMatrix();
    }

    void ToStream(std::ostream &os) { os << *this; }
};

std::istream &operator>>(std::istream &is, std::unique_ptr<SingleOmegaModel> &m) {
    std::string model_name;
    std::string datafile;
    std::string treefile;
    int blmode, nucmode;

    is >> model_name;
    if (model_name != "SingleOmega") {
        std::cerr << "Expected SingleOmega for model name, got " << model_name << "\n";
        exit(1);
    }
    is >> datafile;
    is >> treefile;
    is >> blmode >> nucmode;
    m = std::make_unique<SingleOmegaModel>(datafile, treefile);
    Tracer tracer{*m};
    tracer.read_line(is);
    m->Update();
    return is;
}

std::ostream &operator<<(std::ostream &os, SingleOmegaModel &m) {
    Tracer tracer{m};
    os << "SingleOmega" << '\t';
    os << m.datafile << '\t';
    os << m.treefile << '\t';
    os << m.blmode << '\t';
    os << m.nucmode << '\t';
    tracer.write_line(os);
    return os;
}
