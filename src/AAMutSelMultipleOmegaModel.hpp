#include "AAMutSelCodonMatrixBidimArray.hpp"
#include "Chrono.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "Move.hpp"
#include "MultinomialAllocationVector.hpp"
#include "Permutation.hpp"
#include "PhyloProcess.hpp"
#include "StickBreakingProcess.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"
#include "components/model_decl_utils.hpp"
#include "lib/Dirichlet.hpp"

/**
 * \brief The mutation-selection model with constant fitness landscape over the
 * tree -- double Dirichlet process version.
 *
 * The model is parameterized by
 * - an fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for
 * j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process
 * (assumed homogeneous across sites and lineages)
 * - an array of site-specific amino-acid fitness profiles F_ia, for site i and
 * amino-acid a
 * - an omega multiplier, capturing deviations of the non-syn rate from the
 * model (see Rodrigue and Lartillot, 2107); this parameter is fixed to 1 by
 * default.
 *
 * Site-specific amino-acid fitness profiles are drawn from a Dirichlet process,
 * implemented using a truncated stick-breaking process, of concentration
 * parameter kappa, and truncated at Ncat.
 *
 * The base distribution of this Dirichlet process, G_0, is itself a (truncated
 * stick breaking) mixture of Dirichlet distributions, parameterized by
 * basekappa, and truncated at baseNcat. Each component of this mixture
 * Dirichlet is parameterized by a center (a 20-freqeuncy vector) and a
 * concentration.
 *
 * Concerning the first-level stick-breaking process, by default, Ncat == 100
 * (can be changed using the  -ncat option). As for baseNcat, it is equal to 1,
 * in which case the base distribution G_0 is just a simple Dirichlet (as in
 * Rodrigue et al, 2010). This simple model is probably the one that should be
 * used by default for single-gene analyses. The more complex model with
 * baseNcat > 1 is meant to be used in a multi-gene context (although, even in
 * that case, mixing is still challenging, and not sure whether this model
 * brings important improvement in the end). Thus, baseNcat = 1 appears to be
 * the most reasonable model settings for now.
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 * - omega: fixed to 1 or exponential of rate 1
 * - kappa: exponential of rate 0.1
 * - center of base distribution: uniform Dirichlet
 * - concentration of base distribution: exponential of mean 20.
 *
 * In a multi-gene context, shrinkage across genes can be applied to branch
 * lengths, omega, nucleotide rate parameters (rho and pi), and to the
 * parameters of the base distribution (center and concentration) -- see
 * MultiGeneAAMutSelDSBDPModel.
 *
 */

class AAMutSelMultipleOmegaModel : public ChainComponent {
    std::string datafile, treefile;

    std::unique_ptr<Tracer> tracer;
    std::unique_ptr<const Tree> tree;

    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    // Branch lengths
    double blhypermean;
    double blhyperinvshape;
    SimpleBranchArray<double> *blhypermeanarray;
    GammaWhiteNoise *branchlength;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    // nucleotide rates hyperparameters
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucstat;
    std::vector<double> nucrelrate;
    GTRSubMatrix *nucmatrix;

    int omegaNcat;
    MultinomialAllocationVector *omega_alloc;
    double omega_weight_kappa;
    Dirichlet *omega_weight;

    double omega_shift;
    double delta_omegahypermean;
    double delta_omegahyperinvshape;
    // delta_omega_array of mean delta_omegahypermean and inverse shape parameter
    // delta_omegahyperinvshape
    IIDGamma *delta_omega_array;

    OmegaPathSuffStatArray *siteomegapathsuffstatarray;
    OmegaPathSuffStatArray *componentomegapathsuffstatarray;

    // base distribution G0 is itself a stick-breaking mixture of Dirichlet
    // distributions
    int baseNcat;
    int basemin;
    double basekappa;
    StickBreakingProcess *baseweight;

    OccupancySuffStat *baseprofile_occupancy;
    std::vector<double> basecenterhypercenter;
    double basecenterhyperinvconc;

    IIDDirichlet *basecenterarray;
    double baseconchypermean;
    double baseconchyperinvshape;

    IIDGamma *baseconcentrationarray;
    MultinomialAllocationVector *componentalloc;
    MixtureSelector<std::vector<double>> *componentcenterarray;
    MixtureSelector<double> *componentconcentrationarray;

    // aa fitness arrays across sites are a SBDP process of base G0 defined above
    int Ncat;
    double kappa;
    StickBreakingProcess *weight;

    OccupancySuffStat *profile_occupancy;
    OccupancySuffStat *omega_occupancy;

    MultiDirichlet *componentaafitnessarray;
    DirichletSuffStatArray *basesuffstatarray;

    MultinomialAllocationVector *profile_alloc;

    // an array of codon matrices (one for each distinct aa fitness profile)
    AAMutSelCodonMatrixBidimArray *componentcodonmatrixbidimarray;

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
    DoubleMixtureSelector<SubMatrix> *sitesubmatrixarray;
    DoubleMixtureSelector<AAMutSelOmegaCodonSubMatrix> *sitecodonsubmatrixarray;

    MixtureSelector<std::vector<double>> *siteaafitnessarray;

    PhyloProcess *phyloprocess;

    PathSuffStatArray *sitepathsuffstatarray;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    // currently: shared across genes
    int blmode;
    // currently, free without shrinkage: shared across genes
    int nucmode;
    // currently, shared across genes.
    // free without shrinkage, only with baseNcat = 1
    // free with shrinkage: not really interesting
    int basemode;

    // 0: simple gamma prior
    // 1: mix of (1-pi) at 1 and pi at 1+d, with d ~Gamma(dposomhypermean,dposomhyperinvshape)
    // 2: mix of 2 (modal) gamma distributions: one at 1 and another one with mean > 1
    int omegamode;

    bool flatfitness;

    Chrono aachrono;
    Chrono basechrono;
    Chrono totchrono;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions,
    //! such as specified by branch names)
    //! - inomegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    //! - Ncat: truncation of the first-level stick-breaking process (by default:
    //! 100)
    //! - baseNcat: truncation of the second-level stick-breaking process (by
    //! default: 1)
    AAMutSelMultipleOmegaModel(std::string indatafile, std::string intreefile, int inomegamode,
        int inNcat, int inbaseNcat, int inomegaNcat, double inomegashift, bool inflatfitness)
        : datafile(indatafile),
          treefile(intreefile),
          omegaNcat(inomegaNcat),
          omega_shift(inomegashift),
          baseNcat(inbaseNcat),
          Ncat(inNcat),
          omegamode(inomegamode),
          flatfitness(inflatfitness) {
        init();
        Update();
    }

    virtual ~AAMutSelMultipleOmegaModel() = default;

    void init() {
        blmode = 0;
        nucmode = 0;
        basemode = 0;

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        if (Ncat == -1) { Ncat = Nsite; }
        if (Ncat > Nsite) { Ncat = Nsite; }
        if (Ncat > 100) { Ncat = 100; }

        basemin = 0;
        if (baseNcat < 0) {
            basemin = 1;
            baseNcat = -baseNcat;
            if (baseNcat != 2) {
                std::cerr << "error in basencat\n";
                exit(1);
            }
        }

        if (flatfitness) {
            baseNcat = 1;
            Ncat = 1;
        }

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        std::cerr << "Ncat : " << Ncat << '\n';
        std::cerr << "baseNcat : " << baseNcat << '\n';
        std::cerr << "omegaNcat : " << omegaNcat << '\n';
        std::cerr << "OmegaShift : " << omega_shift << '\n';

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        Nbranch = tree->nb_branches();

        Allocate();
        tracer = std::unique_ptr<Tracer>(new Tracer(*this));
    }

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        if (blmode < 2) { model_node(info, "branchlength", *branchlength); }
        if (nucmode < 2) {
            model_node(info, "nucrelrate", nucrelrate);
            model_node(info, "nucstat", nucstat);
        }
        if (basemode < 2) {
            model_node(info, "basekappa", basekappa);
            model_node(info, "baseweight", *baseweight);
            model_node(info, "componentalloc", *componentalloc);
            model_node(info, "*basecenterarray", *basecenterarray);
            model_node(info, "*baseconcentrationarray", *baseconcentrationarray);
        }
        model_node(info, "kappa", kappa);
        model_node(info, "weight", *weight);
        model_node(info, "componentaafitnessarray", *componentaafitnessarray);
        model_node(info, "profile_alloc", *profile_alloc);
        model_node(info, "delta_omega_array", *delta_omega_array);
        model_node(info, "omega_alloc", *omega_alloc);
        model_node(info, "omega_weight", *omega_weight);
        model_node(info, "delta_omegahyperinvshape", delta_omegahyperinvshape);
        model_node(info, "delta_omegahypermean", delta_omegahypermean);


        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        // 3x: per coding site (and not per nucleotide site)
        model_stat(info, "length", [this]() { return 3 * branchlength->GetTotalLength(); });
        model_stat(info, "dnds", [this]() { return GetPredictedEffectivedNdS(); });
        model_stat(info, "omega0", [this]() { return GetPredictedOmegaKnot(); });
        model_stat(info, "omega", [this]() { return GetMeanOmega(); });
        model_stat(
            info, "omegaent", [this]() { return Random::GetEntropy(omega_weight->GetArray()); });
        model_stat(info, "ncluster", [this]() { return GetNcluster(); });
        model_stat(info, "kappa", kappa);
        if (baseNcat > 1) {
            if (basemin) {
                model_stat(info, "basencluster", [this]() { return GetBaseNcluster(); });
                model_stat(info, "baseweight1", [this]() { return baseweight->GetVal(1); });
            } else {
                model_stat(info, "basencluster", [this]() { return GetBaseNcluster(); });
                model_stat(info, "basekappa", basekappa);
            }
        }
        model_stat(info, "aaent", [this]() { return GetMeanAAEntropy(); });
        model_stat(info, "meanaaconc", [this]() { return GetMeanComponentAAConcentration(); });
        model_stat(info, "aacenterent", [this]() { return GetMeanComponentAAEntropy(); });
        model_stat(info, "statent", [this]() { return GetNucRREntropy(); });
        model_stat(info, "rrent", [this]() { return GetNucStatEntropy(); });
    }

    //! allocate the model (data structures)
    void Allocate() {
        // Branch lengths
        blhypermean = 0.1;
        blhyperinvshape = 1.0;
        blhypermeanarray = new SimpleBranchArray<double>(*tree, blhypermean);
        branchlength = new GammaWhiteNoise(*tree, *blhypermeanarray, blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, std::vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, std::vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));
        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        // base distribution (can be skipped)
        basekappa = 1.0;
        baseweight = new StickBreakingProcess(baseNcat, basekappa);
        baseprofile_occupancy = new OccupancySuffStat(baseNcat);

        basecenterhypercenter.assign(Naa, 1.0 / Naa);
        basecenterhyperinvconc = 1.0 / Naa;

        basecenterarray = new IIDDirichlet(baseNcat, basecenterhypercenter, basecenterhyperinvconc);
        basecenterarray->SetUniform();

        baseconchypermean = Naa;
        baseconchyperinvshape = 1.0;
        baseconcentrationarray = new IIDGamma(baseNcat, baseconchypermean, baseconchyperinvshape);
        for (int k = 0; k < baseNcat; k++) { (*baseconcentrationarray)[k] = 20.0; }
        if (basemin == 1) { (*baseconcentrationarray)[0] = 1.0; }
        // suff stats for component aa fitness arrays
        basesuffstatarray = new DirichletSuffStatArray(baseNcat, Naa);
        componentalloc = new MultinomialAllocationVector(Ncat, baseweight->GetArray());
        componentcenterarray =
            new MixtureSelector<std::vector<double>>(basecenterarray, componentalloc);
        componentconcentrationarray =
            new MixtureSelector<double>(baseconcentrationarray, componentalloc);

        //
        // (truncated) Dirichlet mixture of aa fitness profiles
        //

        // Ncat fitness profiles iid from the base distribution
        componentaafitnessarray =
            new MultiDirichlet(componentcenterarray, componentconcentrationarray);
        if (flatfitness) { componentaafitnessarray->Flatten(); }

        // mixture weights (truncated stick breaking process)
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat, kappa);

        // site allocations to the mixture (multinomial allocation)
        profile_alloc = new MultinomialAllocationVector(Nsite, weight->GetArray());

        // occupancy suff stats of site allocations for resampling weights (for both profile and
        // omega)
        profile_occupancy = new OccupancySuffStat(Ncat);
        omega_occupancy = new OccupancySuffStat(omegaNcat);

        // omega (fixed to 1 by default)
        delta_omegahypermean = 1.0;
        delta_omegahyperinvshape = 1.0;
        delta_omega_array = new IIDGamma(omegaNcat, delta_omegahypermean, delta_omegahyperinvshape);

        // will be a constant
        omega_weight_kappa = 1.0;
        omega_weight = new Dirichlet(omegaNcat, omega_weight_kappa);
        omega_alloc = new MultinomialAllocationVector(GetNsite(), omega_weight->GetArray());

        // Ncat*omegaNcat mut sel codon matrices (based on the Ncat fitness profiles of the mixture,
        // and the omegaNcat finite mixture for omega)
        componentcodonmatrixbidimarray = new AAMutSelCodonMatrixBidimArray(GetCodonStateSpace(),
            nucmatrix, componentaafitnessarray, delta_omega_array, omega_shift);

        // Bidimselector, specifying which codon matrix should be used for each site
        sitesubmatrixarray = new DoubleMixtureSelector<SubMatrix>(
            componentcodonmatrixbidimarray, profile_alloc, omega_alloc);

        // Bidimselector, specifying which codon matrix should be used for each site (needed to
        // collect omegapathsuffstatarray)
        sitecodonsubmatrixarray = new DoubleMixtureSelector<AAMutSelOmegaCodonSubMatrix>(
            componentcodonmatrixbidimarray, profile_alloc, omega_alloc);

        // selector, specifying which aa fitness array should be used for each site
        siteaafitnessarray =
            new MixtureSelector<std::vector<double>>(componentaafitnessarray, profile_alloc);

        phyloprocess =
            new PhyloProcess(tree.get(), codondata, branchlength, nullptr, sitesubmatrixarray);
        phyloprocess->Unfold();

        sitepathsuffstatarray = new PathSuffStatArray(Nsite);

        siteomegapathsuffstatarray = new OmegaPathSuffStatArray(Nsite);
        componentomegapathsuffstatarray = new OmegaPathSuffStatArray(omegaNcat);
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //! return current omega value for omega mixture of component k
    double GetComponentOmega(int k) const { return omega_shift + delta_omega_array->GetVal(k); }

    //! return current omega value for omega mixture of site
    double GetSiteOmega(int site) const { return GetComponentOmega(omega_alloc->GetVal(site)); }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void CorruptCodonMatrices() { componentcodonmatrixbidimarray->CorruptCodonMatrices(); }

    //! \brief tell codon matrices (of profile component k) that its parameters have changed and
    //! that it should be updated
    void CorruptProfileCodonMatrices(int k) {
        componentcodonmatrixbidimarray->CorruptRowCodonMatrices(k);
    }

    //! \brief tell codon matrices (of omega component k) that its parameters have changed and
    //! that it should be updated
    void CorruptOmegaCodonMatrix(int k) {
        componentcodonmatrixbidimarray->CorruptRowCodonMatrices(k);
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
        UpdateNucMatrix();
        CorruptCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void Update() {
        UpdateBaseOccupancies();
        UpdateProfileOccupancies();
        UpdateOmegaOccupancies();
        UpdateMatrices();
        ResampleSub(1.0);
    }

    void PostPred(std::string name) {
        UpdateBaseOccupancies();
        UpdateOmegaOccupancies();
        UpdateProfileOccupancies();
        UpdateMatrices();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior (up to some constant)
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) { total += BranchLengthsLogPrior(); }
        if (nucmode < 2) { total += NucRatesLogPrior(); }
        if (basemode < 2) {
            if (baseNcat > 1) {
                total += BaseStickBreakingHyperLogPrior();
                total += BaseStickBreakingLogPrior();
            }
            total += BaseLogPrior();
        }
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AALogPrior();
        if (omegamode < 2) { total += DeltaOmegaLogProb(); }
        return total;
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const {
        double ret = branchlength->GetLogProb();
        return ret;
    }

    //! \brief log prior over delta_omega (for all components) + HyperLogPrior
    double DeltaOmegaLogProb() const { return DeltaOmegaLogPrior() + DeltaOmegaHyperLogPrior(); }

    //! \brief log prior over delta_omega (for all components)
    double DeltaOmegaLogPrior() const { return delta_omega_array->GetLogProb(); }

    //! \brief log prior over delta_omega (for component i)
    double DeltaOmegaLogPrior(int i) const { return delta_omega_array->GetLogProb(i); }

    //! \brief log prior over hyperparameter of prior over omega (here,
    //! delta_omegahyperinvshape ~ exponential of rate 10)
    //! delta_omegahypermean ~ exponential of rate 10)
    double DeltaOmegaHyperLogPrior() const {
        return -(delta_omegahyperinvshape + delta_omegahypermean) / 10;
    }

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(
            nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //! log prior over concentration parameters basekappa of mixture of base
    //! distribution
    double BaseStickBreakingHyperLogPrior() const { return -basekappa / 10; }

    //! log prior over weights of stick breaking process of base distribution
    double BaseStickBreakingLogPrior() const { return baseweight->GetLogProb(basekappa); }

    //! log prior over concentration parameters kappa of stick-breaking mixture of
    //! amino-acid profiles
    double StickBreakingHyperLogPrior() const { return -kappa / 10; }

    //! log prior over weights of stick breaking process of amino-acid profiles
    double StickBreakingLogPrior() const { return weight->GetLogProb(kappa); }

    //! log prior over base center and concentration parameters
    double BaseLogPrior() const {
        double total = 0;
        total += basecenterarray->GetLogProb();
        total += baseconcentrationarray->GetLogProb();
        if (std::isinf(total)) {
            std::cerr << "in BaseLogPrior: inf\n";
            exit(1);
        }
        return total;
    }

    //! log prior over base center and concentration parameters of component k of
    //! base distribution
    double BaseLogPrior(int k) const {
        double total = 0;
        total += basecenterarray->GetLogProb(k);
        total += baseconcentrationarray->GetLogProb(k);
        return total;
    }

    //! log prior of amino-acid fitness profiles
    double AALogPrior() const { return componentaafitnessarray->GetLogProb(); }

    //! log prior of amino-acid fitness profile k
    double AALogPrior(int k) const { return componentaafitnessarray->GetLogProb(k); }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        double tot{0};
        for (int site{0}; site < GetNsite(); site++) {
            tot += sitepathsuffstatarray->GetVal(site).GetLogProb(sitesubmatrixarray->GetVal(site));
        }
        return tot;
    }

    //! return log prob of the substitution mappings over sites allocated to
    //! component k of the mixture
    double PathSuffStatLogProb(int k) const {
        double tot{0};
        for (int site{0}; site < GetNsite(); site++) {
            if (profile_alloc->GetVal(site) == k) {
                tot += sitepathsuffstatarray->GetVal(site).GetLogProb(
                    sitesubmatrixarray->GetVal(site));
            }
        }
        return tot;
    }

    //! return log prob of the substitution mappings over sites allocated to omega
    //! component k of the omega mixture
    double PathSuffStatOmegaLogProb(int k) const {
        return componentomegapathsuffstatarray->GetVal(k).GetLogProb(GetComponentOmega(k));
    }

    //! return log prob of first-level mixture components (i.e. all amino-acid
    //! profiles drawn from component k of the base distribution), as a function
    //! of the center and concentration parameters of this component
    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(
            basecenterarray->GetVal(k), 1.0 / baseconcentrationarray->GetVal(k));
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! log prob factor to be recomputed when moving nucleotide mutation rate
    //! parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + PathSuffStatLogProb(); }

    //! log prob factor to be recomputed when moving aa hyper params (center and
    //! concentration) for component k of base distribution
    double BaseLogProb(int k) const { return BaseLogPrior(k) + BaseSuffStatLogProb(k); }

    //! log prob factor to be recomputed when moving basekappa, the concentration
    //! parameter of the second-level strick breaking process (base distribution)
    double BaseStickBreakingHyperLogProb() const {
        return BaseStickBreakingHyperLogPrior() + BaseStickBreakingLogPrior();
    }

    //! log prob factor to be recomputed when moving kappa, the concentration
    //! parameter of the first-level strick breaking process
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics of substitution mappings across sites
    void CollectSitePathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    //! collect omega sufficient statistics of substitution mappings across sites
    void CollectSiteOmegaPathSuffStat() {
        siteomegapathsuffstatarray->Clear();
        siteomegapathsuffstatarray->AddSuffStat(*sitecodonsubmatrixarray, *sitepathsuffstatarray);
    }

    //! gather site-specific omega sufficient statistics component-wise (per omega component)
    void CollectComponentOmegaPathSuffStat() {
        componentomegapathsuffstatarray->Clear();
        componentomegapathsuffstatarray->Add(*siteomegapathsuffstatarray, *omega_alloc);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //-------------------
    //  Moves
    //-------------------

    //! complete MCMC move schedule
    double Move() {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        UpdateMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            totchrono.Start();
            if (blmode < 2) { MoveBranchLengths(); }

            CollectSitePathSuffStat();

            if (nucmode < 2) { MoveNucRates(); }

            if (omegamode < 2) { MoveOmegaMixture(3); }

            if (!flatfitness) {
                aachrono.Start();
                MoveAAMixture(3);
                aachrono.Stop();

                basechrono.Start();
                if (basemode < 2) { MoveBase(3); }
                basechrono.Stop();
            }

            totchrono.Stop();
        }
    }

    //! MH move on base mixture
    void MoveBase(int nrep) {
        if (baseNcat > 1) { ResampleBaseAlloc(); }
        MoveBaseMixture(nrep);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MCMC move schedule on branch lengths
    void MoveBranchLengths() { ResampleBranchLengths(); }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        Move::Profile(nucrelrate, 0.1, 1, 3, &AAMutSelMultipleOmegaModel::NucRatesLogProb,
            &AAMutSelMultipleOmegaModel::UpdateMatrices, this);
        Move::Profile(nucrelrate, 0.03, 3, 3, &AAMutSelMultipleOmegaModel::NucRatesLogProb,
            &AAMutSelMultipleOmegaModel::UpdateMatrices, this);
        Move::Profile(nucrelrate, 0.01, 3, 3, &AAMutSelMultipleOmegaModel::NucRatesLogProb,
            &AAMutSelMultipleOmegaModel::UpdateMatrices, this);

        Move::Profile(nucstat, 0.1, 1, 3, &AAMutSelMultipleOmegaModel::NucRatesLogProb,
            &AAMutSelMultipleOmegaModel::UpdateMatrices, this);
        Move::Profile(nucstat, 0.01, 1, 3, &AAMutSelMultipleOmegaModel::NucRatesLogProb,
            &AAMutSelMultipleOmegaModel::UpdateMatrices, this);
    }

    //! MCMC module for the mixture amino-acid fitness profiles
    void MoveAAMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveAAProfiles();
            ResampleEmptyProfileComponents();
            ResampleProfileAlloc();
            ProfileLabelSwitchingMove();
            ResampleProfileWeights();
            MoveKappa();
            CorruptCodonMatrices();
        }
    }

    //! resample empty components of the mixture from prior
    void ResampleEmptyProfileComponents() {
        componentaafitnessarray->PriorResample(*profile_occupancy);
        componentcodonmatrixbidimarray->CorruptCodonMatricesRowOccupancy(*profile_occupancy);
    }

    //! MH move on amino-acid fitness profiles (occupied components only)
    void MoveAAProfiles() {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double CompMoveAAProfiles(int nrep) {
        MoveAA(1.0, 1, nrep);
        MoveAA(0.1, 3, nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MulMoveAAProfiles(int nrep) {
        MoveAAGamma(3.0, nrep);
        MoveAAGamma(1.0, nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double MoveAA(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Naa];
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) {
                std::vector<double> &aa = (*componentaafitnessarray)[i];
                for (int rep = 0; rep < nrep; rep++) {
                    for (int l = 0; l < Naa; l++) { bk[l] = aa[l]; }
                    double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                    deltalogprob += loghastings;
                    CorruptProfileCodonMatrices(i);
                    deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        for (int l = 0; l < Naa; l++) { aa[l] = bk[l]; }
                        CorruptProfileCodonMatrices(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! helper function: log density of 20 gammas
    double GammaAALogPrior(
        const std::vector<double> &x, const std::vector<double> &aacenter, double aaconc) {
        double total = 0;
        for (int l = 0; l < Naa; l++) {
            total += (aaconc * aacenter[l] - 1) * log(x[l]) - x[l] -
                     Random::logGamma(aaconc * aacenter[l]);
        }
        return total;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MoveAAGamma(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) {
                double aaconc = componentconcentrationarray->GetVal(i);
                const std::vector<double> &aacenter = componentcenterarray->GetVal(i);

                std::vector<double> &aa = (*componentaafitnessarray)[i];
                std::vector<double> x(Naa, 0);
                double z = Random::sGamma(aaconc);
                for (int l = 0; l < Naa; l++) { x[l] = z * aa[l]; }

                double bkz = z;
                std::vector<double> bkx = x;
                std::vector<double> bkaa = aa;

                for (int rep = 0; rep < nrep; rep++) {
                    double deltalogprob =
                        -GammaAALogPrior(x, aacenter, aaconc) - PathSuffStatLogProb(i);

                    double loghastings = 0;
                    z = 0;
                    for (int l = 0; l < Naa; l++) {
                        double m = tuning * (Random::Uniform() - 0.5);
                        double e = exp(m);
                        x[l] *= e;
                        z += x[l];
                        loghastings += m;
                    }
                    for (int l = 0; l < Naa; l++) {
                        aa[l] = x[l] / z;
                        if (aa[l] < 1e-50) { aa[l] = 1e-50; }
                    }

                    deltalogprob += loghastings;

                    CorruptProfileCodonMatrices(i);

                    deltalogprob += GammaAALogPrior(x, aacenter, aaconc) + PathSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                        bkaa = aa;
                        bkx = x;
                        bkz = z;
                    } else {
                        aa = bkaa;
                        x = bkx;
                        z = bkz;
                        CorruptProfileCodonMatrices(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! Gibbs resample mixture allocations
    void ResampleProfileAlloc() {
        std::vector<double> postprob(Ncat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetProfileAllocPostProb(i, postprob);
            profile_alloc->GibbsResample(i, postprob);
        }
        UpdateProfileOccupancies();
    }

    //! update mixture profile occupancy suff stats (for resampling mixture weights)
    void UpdateProfileOccupancies() {
        profile_occupancy->Clear();
        profile_occupancy->AddSuffStat(*profile_alloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetProfileAllocPostProb(int site, std::vector<double> &postprob) {
        double max = 0;
        const std::vector<double> &w = weight->GetArray();
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        for (int i = 0; i < Ncat; i++) {
            double tmp = suffstat.GetLogProb(
                componentcodonmatrixbidimarray->GetVal(i, omega_alloc->GetVal(site)));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) { max = tmp; }
        }

        double total = 0;
        for (int i = 0; i < Ncat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < Ncat; i++) { postprob[i] /= total; }
    }

    //! MCMC sequence for label switching moves
    void ProfileLabelSwitchingMove() {
        Permutation permut(Ncat);
        weight->LabelSwitchingMove(5, *profile_occupancy, permut);
        profile_alloc->Permute(permut);
        componentaafitnessarray->Permute(permut);
    }

    //! Gibbs resample mixture weights (based on profile_occupancy suff stats)
    void ResampleProfileWeights() { weight->GibbsResample(*profile_occupancy); }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa() {
        Move::Scaling(kappa, 1.0, 10, &AAMutSelMultipleOmegaModel::StickBreakingHyperLogProb,
            &AAMutSelMultipleOmegaModel::NoUpdate, this);
        Move::Scaling(kappa, 0.3, 10, &AAMutSelMultipleOmegaModel::StickBreakingHyperLogProb,
            &AAMutSelMultipleOmegaModel::NoUpdate, this);
    }

    //! MCMC module for the base mixture
    void MoveBaseMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseComponents(10);
            ResampleBaseEmptyComponents();
            if (baseNcat > 1) {
                if (!basemin) { BaseLabelSwitchingMove(); }
                ResampleBaseWeights();
                MoveBaseKappa();
            }
        }
    }

    //! MCMC module for moving the center and concentration parameters of the
    //! components of the the base mixture
    void MoveBaseComponents(int nrep) {
        CollectBaseSuffStat();
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseCenters(1.0, 1);
            MoveBaseCenters(1.0, 3);
            MoveBaseCenters(0.3, 3);
            MoveBaseConcentrations(1.0);
            MoveBaseConcentrations(0.3);
        }
    }

    //! collect suff stats for moving center and concentration parameters of the
    //! base mixture
    void CollectBaseSuffStat() {
        basesuffstatarray->Clear();
        componentaafitnessarray->AddSuffStat(*basesuffstatarray, *componentalloc);
    }

    //! MCMC module for moving the center parameters of the components of the the
    //! base mixture
    double MoveBaseCenters(double tuning, int n) {
        double nacc = 0;
        double ntot = 0;
        std::vector<double> bk(Naa, 0);
        for (int k = 0; k < baseNcat; k++) {
            if (baseprofile_occupancy->GetVal(k)) {
                std::vector<double> &aa = (*basecenterarray)[k];
                bk = aa;
                double deltalogprob = -BaseLogProb(k);
                double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                deltalogprob += loghastings;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    aa = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MCMC module for moving the concentration parameters of the components of
    //! the the base mixture
    double MoveBaseConcentrations(double tuning) {
        double nacc = 0;
        double ntot = 0;
        for (int k = basemin; k < baseNcat; k++) {
            if (baseprofile_occupancy->GetVal(k)) {
                double &c = (*baseconcentrationarray)[k];
                double bk = c;
                double deltalogprob = -BaseLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                c *= e;
                deltalogprob += m;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    c = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! resample empty components of the base mixture from the prior
    void ResampleBaseEmptyComponents() {
        basecenterarray->PriorResample(*baseprofile_occupancy);
        baseconcentrationarray->PriorResample(*baseprofile_occupancy);
        if (basemin == 1) { (*baseconcentrationarray)[0] = 1.0; }
    }

    //! Gibbs resample base mixture allocations
    void ResampleBaseAlloc() {
        std::vector<double> postprob(baseNcat, 0);
        for (int i = 0; i < Ncat; i++) {
            GetBaseAllocPostProb(i, postprob);
            componentalloc->GibbsResample(i, postprob);
            if ((componentalloc->GetVal(i) < 0) || (componentalloc->GetVal(i) >= baseNcat)) {
                std::cerr << "error in ResampleBaseAlloc: out of bound\n";
                exit(1);
            }
        }
        UpdateBaseOccupancies();
    }

    //! update base profile_occupancy suff stats (for moving base weights)
    void UpdateBaseOccupancies() {
        baseprofile_occupancy->Clear();
        baseprofile_occupancy->AddSuffStat(*componentalloc);
    }

    //! get allocation posterior probability of a component of the first-level
    //! mixture to the components of the second-level mixture
    void GetBaseAllocPostProb(int cat, std::vector<double> &postprob) {
        double max = 0;
        const std::vector<double> &w = baseweight->GetArray();
        for (int i = 0; i < baseNcat; i++) {
            double tmp = Random::logDirichletDensity(componentaafitnessarray->GetVal(cat),
                basecenterarray->GetVal(i), baseconcentrationarray->GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) { max = tmp; }
        }

        double total = 0;
        for (int i = 0; i < baseNcat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < baseNcat; i++) { postprob[i] /= total; }
    }

    //! MCMC sequence for label switching moves of the base mixture
    void BaseLabelSwitchingMove() {
        Permutation permut(baseNcat);
        baseweight->LabelSwitchingMove(5, *baseprofile_occupancy, permut);
        componentalloc->Permute(permut);
        basecenterarray->Permute(permut);
        baseconcentrationarray->Permute(permut);
        basesuffstatarray->Permute(permut);
    }

    //! Gibbs resample base mixture weights (based on profile_occupancy suff stats)
    void ResampleBaseWeights() { baseweight->GibbsResample(*baseprofile_occupancy); }

    //! MH move on basekappa, concentration parameter of the base mixture
    void MoveBaseKappa() {
        Move::Scaling(basekappa, 1.0, 10,
            &AAMutSelMultipleOmegaModel::BaseStickBreakingHyperLogProb,
            &AAMutSelMultipleOmegaModel::NoUpdate, this);
        Move::Scaling(basekappa, 0.3, 10,
            &AAMutSelMultipleOmegaModel::BaseStickBreakingHyperLogProb,
            &AAMutSelMultipleOmegaModel::NoUpdate, this);
    }

    //! MCMC module for the mixture of omega
    void MoveOmegaMixture(int nrep) {
        CollectSiteOmegaPathSuffStat();
        for (int rep = 0; rep < nrep; rep++) {
            ResampleOmegaAlloc();
            UpdateOmegaOccupancies();
            CollectComponentOmegaPathSuffStat();
            MoveOmegaValues(1.0, 3);
            MoveOmegaValues(0.1, 3);
            ResampleEmptyOmegaComponents();
            ResampleOmegaWeights();
            MoveOmegaHyper();
        }
        CorruptCodonMatrices();
    }

    //! MH move on omega (per component): scaling move
    double MoveOmegaValues(double tuning, int nrep) {
        double bk{0}, nacc{0}, ntot{0};
        for (int i = 0; i < omegaNcat; i++) {
            if (omega_occupancy->GetVal(i)) {
                for (int rep = 0; rep < nrep; rep++) {
                    bk = delta_omega_array->GetVal(i);

                    double deltalogprob = -(PathSuffStatOmegaLogProb(i) + DeltaOmegaLogPrior(i));
                    double loghastings = tuning * (Random::Uniform() - 0.5);
                    double hastings = exp(loghastings);

                    (*delta_omega_array)[i] *= hastings;

                    deltalogprob += loghastings;
                    deltalogprob += PathSuffStatOmegaLogProb(i) + DeltaOmegaLogPrior(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        (*delta_omega_array)[i] = bk;
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    void ResampleEmptyOmegaComponents() { delta_omega_array->PriorResample(*omega_occupancy); };

    //! update mixture omega occupancy suff stats (for resampling mixture weights)
    void UpdateOmegaOccupancies() {
        omega_occupancy->Clear();
        omega_occupancy->AddSuffStat(*omega_alloc);
    }

    //! Gibbs resample omega mixture allocations
    void ResampleOmegaAlloc() {
        std::vector<double> postprob(omegaNcat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetOmegaAllocPostProb(i, postprob);
            omega_alloc->GibbsResample(i, postprob);
        }
    }

    //! get omega allocation posterior probabilities for a given site
    double GetOmegaAllocPostProb(int site, std::vector<double> &postprob) {
        double max = 0;

        const std::vector<double> &w = omega_weight->GetArray();
        const OmegaPathSuffStat &suffstat = siteomegapathsuffstatarray->GetVal(site);

        for (int i = 0; i < omegaNcat; i++) {
            double tmp = suffstat.GetLogProb(GetComponentOmega(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) { max = tmp; }
        }

        double total = 0;
        for (int i = 0; i < omegaNcat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < omegaNcat; i++) { postprob[i] /= total; }

        return log(total) + max;
    }

    //! Gibbs resample omega mixture weights (based on omega occupancy suff stats)
    void ResampleOmegaWeights() { omega_weight->GibbsResample(*omega_occupancy); }


    //! Scaling moves on omega hyperparameters
    void MoveOmegaHyper() {
        if (omegaNcat > 1) {
            Move::Scaling(delta_omegahypermean, 1.0, 10,
                &AAMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &AAMutSelMultipleOmegaModel::NoUpdate, this);
            Move::Scaling(delta_omegahypermean, 0.3, 10,
                &AAMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &AAMutSelMultipleOmegaModel::NoUpdate, this);
            Move::Scaling(delta_omegahyperinvshape, 1.0, 10,
                &AAMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &AAMutSelMultipleOmegaModel::NoUpdate, this);
            Move::Scaling(delta_omegahyperinvshape, 0.3, 10,
                &AAMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &AAMutSelMultipleOmegaModel::NoUpdate, this);
        }
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    //! return number of occupied components in first-level mixture (mixture of
    //! amino-acid fitness profiles)
    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) { n++; }
        }
        return n;
    }

    //! return number of occupied components in base distribution
    int GetBaseNcluster() const {
        int n = 0;
        for (int i = 0; i < baseNcat; i++) {
            if (baseprofile_occupancy->GetVal(i)) { n++; }
        }
        return n;
    }


    //! return mean omega
    double GetMeanOmega() const {
        double tot = 0;
        for (int i{0}; i < omegaNcat; i++) { tot += GetComponentOmega(i); }
        return tot / omegaNcat;
    }

    //! return mean entropy of amino-acd fitness profiles
    double GetMeanAAEntropy() const { return componentaafitnessarray->GetMeanEntropy(); }

    //! return mean of concentration parameters of base distribution
    double GetMeanComponentAAConcentration() const {
        double tot = 0;
        double totw = 0;
        for (int i = basemin; i < baseNcat; i++) {
            tot += baseprofile_occupancy->GetVal(i) * baseconcentrationarray->GetVal(i);
            totw += baseprofile_occupancy->GetVal(i);
        }
        return tot / totw;
    }

    //! return mean entropy of centers of base distribution
    double GetMeanComponentAAEntropy() const {
        double tot = 0;
        for (int i = 0; i < baseNcat; i++) {
            tot +=
                baseprofile_occupancy->GetVal(i) * Random::GetEntropy(basecenterarray->GetVal(i));
        }
        return tot / Ncat;
    }

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return predicted omega induced by the mutation-selection balance
    double GetPredictedOmegaKnot() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += sitecodonsubmatrixarray->GetVal(i).GetPredictedDNDS();
        }
        mean /= GetNsite();
        return mean;
    }

    //! return effective dnds taking into the mutation-selection balance and omega
    double GetPredictedEffectivedNdS() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += GetSiteOmega(i) * sitecodonsubmatrixarray->GetVal(i).GetPredictedDNDS();
        }
        mean /= GetNsite();
        return mean;
    }

    const std::vector<double> &GetProfile(int i) const { return siteaafitnessarray->GetVal(i); }

    void ToStream(std::ostream &os) const {
        os << "AAMutSelMultipleOmega" << '\t';
        os << datafile << '\t' << treefile << '\t';
        os << omegamode << '\t';
        os << Ncat << '\t' << baseNcat << '\t';
        os << omegaNcat << '\t' << omega_shift << '\t';
        tracer->write_line(os);
    }

    AAMutSelMultipleOmegaModel(std::istream &is) {
        std::string model_name;
        is >> model_name;
        if (model_name != "AAMutSelMultipleOmega") {
            std::cerr << "Expected AAMutSelMultipleOmega for model name, got " << model_name
                      << "\n";
            exit(1);
        }
        is >> datafile >> treefile;
        is >> omegamode;
        is >> Ncat >> baseNcat >> omegaNcat >> omega_shift;

        init();
        tracer->read_line(is);
        Update();
    };
};
