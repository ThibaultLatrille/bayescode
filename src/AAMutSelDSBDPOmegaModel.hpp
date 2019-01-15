#include "AAMutSelOmegaCodonSubMatrix.hpp"
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
#include "PolyProcess.hpp"
#include "PolySuffStat.hpp"
#include "StickBreakingProcess.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"

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
 * - branch lengths iid exponential
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

class AAMutSelDSBDPOmegaModel : public ChainComponent {
    std::string datafile, treefile;

    bool polymorphism_aware;

    std::unique_ptr<Tracer> tracer;
    std::unique_ptr<const Tree> tree;

    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;
    PolyData *polydata{nullptr};

    int Nsite;
    int Ntaxa;
    int Nbranch;

    // Branch lengths
    double blhypermean;
    double blhyperinvshape;
    SimpleBranchArray<double> *blhypermeanarray;
    GammaWhiteNoise *branchlength;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    // nucleotide rates hyperparameters
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucstat;
    std::vector<double> nucrelrate;
    GTRSubMatrix *nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
    double omega;
    OmegaPathSuffStat omegapathsuffstat;

    // base distribution G0 is itself a stick-breaking mixture of Dirichlet
    // distributions

    int baseNcat;
    int basemin;
    double basekappa;
    StickBreakingProcess *baseweight;
    OccupancySuffStat *baseoccupancy;

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
    OccupancySuffStat *occupancy;

    MultiDirichlet *componentaafitnessarray;
    DirichletSuffStatArray *basesuffstatarray;

    MultinomialAllocationVector *sitealloc;

    // an array of codon matrices (one for each distinct aa fitness profile)
    AAMutSelOmegaCodonSubMatrixArray *componentcodonmatrixarray;

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
    MixtureSelector<SubMatrix> *sitesubmatrixarray;

    // this one is used by PolyProcess: has to be a Selector<vector<double>>
    MixtureSelector<std::vector<double>> *siteaafitnessarray;

    PhyloProcess *phyloprocess;

    // global theta (4*Ne*u) used for polymorphism
    double theta_scale;
    double thetamax;
    HomogeneousScaledMutationRate *theta;

    PolyProcess *polyprocess{nullptr};
    PoissonRandomField *poissonrandomfield{nullptr};

    PathSuffStatArray *sitepathsuffstatarray;
    PathSuffStatArray *componentpathsuffstatarray;

    PolySuffStatArray *sitepolysuffstatarray{nullptr};
    PolySuffStatArray *componentpolysuffstatarray{nullptr};

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

    // currently: fixed or free with shrinkage
    int omegaprior;

    // 0: simple gamma prior
    // 1: mix of (1-pi) at 1 and pi at 1+d, with d ~Gamma(dposomhypermean,dposomhyperinvshape)
    // 2: mix of 2 (modal) gamma distributions: one at 1 and another one with mean > 1
    int omegamode;

    double dposompi, dposomhypermean, dposomhyperinvshape;

    Chrono aachrono;
    Chrono basechrono;
    Chrono totchrono;

    double acca1, acca2, acca3, acca4;
    double tota1, tota2, tota3, tota4;
    double accb1, accb2, accb3, accb4;
    double totb1, totb2, totb3, totb4;

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
    //! - polymorphism_aware: boolean to force using polymorphism data
    AAMutSelDSBDPOmegaModel(std::string indatafile, std::string intreefile, int inomegamode,
        int inomegaprior, double indposompi, double indposomhypermean, double indposomhyperinvshape,
        int inNcat, int inbaseNcat, bool polymorphism_aware)
        : datafile(indatafile),
          treefile(intreefile),
          polymorphism_aware(polymorphism_aware),
          baseNcat(inbaseNcat),
          Ncat(inNcat),
          omegaprior(inomegaprior),
          omegamode(inomegamode),
          dposompi(indposompi),
          dposomhypermean(indposomhypermean),
          dposomhyperinvshape(indposomhyperinvshape) {
        init();
        Update();
    }

    virtual ~AAMutSelDSBDPOmegaModel() = default;

    void init() {
        blmode = 0;
        nucmode = 0;
        basemode = 0;

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        if (polymorphism_aware) { polydata = new PolyData(codondata, datafile); }

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

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        std::cerr << "ncat : " << Ncat << '\n';
        std::cerr << "basencat : " << baseNcat << '\n';

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        Nbranch = tree->nb_branches();

        // not sure this was not needed
        // tree->SetIndices();

        acca1 = acca2 = acca3 = acca4 = 0;
        tota1 = tota2 = tota3 = tota4 = 0;
        accb1 = accb2 = accb3 = accb4 = 0;
        totb1 = totb2 = totb3 = totb4 = 0;

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
        model_node(info, "sitealloc", *sitealloc);
        if (omegamode < 2) { model_node(info, "omega; ", omega); }
        if (polyprocess != nullptr) { model_node(info, "theta; ", theta); }


        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        // 3x: per coding site (and not per nucleotide site)
        model_stat(info, "length", [this]() { return 3 * branchlength->GetTotalLength(); });
        model_stat(info, "omega", omega);
        model_stat(info, "theta", theta);
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
        model_stat(info, "statent", [&]() { return Random::GetEntropy(nucstat); });
        model_stat(info, "rrent", [&]() { return Random::GetEntropy(nucrelrate); });
    }

    //! \brief set estimation method for branch lengths
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetBLMode(int mode) { blmode = mode; }

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetNucMode(int mode) { nucmode = mode; }

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 3: fixed to 1
    //! - mode == 2: shared and estimated across genes: currently not implemented
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, either mode 3 and mode 0 can be used -- default
    //! mode is 3.
    void SetOmegaMode(int mode) { omegamode = mode; }

    //! \brief set prior for omega
    //
    //! 0: simple gamma prior
    //! 1: mix of (1-pi) at 1 and pi at 1+d, with d ~
    //! Gamma(dposomhypermean,dposomhyperinvshape)
    void SetOmegaPrior(int inprior) { omegaprior = inprior; }

    //! \brief set estimation method for center and concentration parameters of
    //! base distribution
    //!
    //! - mode == 3: shared across genes and fixed to externally given empirical
    //! values: currently not implemented
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage): currently not implemented
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, either mode 3 and mode 0 can be used -- default
    //! mode is 3.
    void SetBaseMode(int mode) { basemode = mode; }

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
        baseoccupancy = new OccupancySuffStat(baseNcat);

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

        // mixture weights (truncated stick breaking process)
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat, kappa);

        // site allocations to the mixture (multinomial allocation)
        sitealloc = new MultinomialAllocationVector(Nsite, weight->GetArray());

        // occupancy suff stats of site allocations (for resampling weights)
        occupancy = new OccupancySuffStat(Ncat);

        // global omega (fixed to 1 by default)
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        omega = 1.0;

        // Ncat mut sel codon matrices (based on the Ncat fitness profiles of the
        // mixture)
        componentcodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(
            GetCodonStateSpace(), nucmatrix, componentaafitnessarray, omega);

        // selector, specifying which codon matrix should be used for each site
        sitesubmatrixarray = new MixtureSelector<SubMatrix>(componentcodonmatrixarray, sitealloc);

        // selector, specifying which aa fitness array should be used for each site
        siteaafitnessarray =
            new MixtureSelector<std::vector<double>>(componentaafitnessarray, sitealloc);

        // global theta (4*Ne*u = 1e-5 by default, and maximum value 0.1)
        theta_scale = 1e-5;
        theta = new HomogeneousScaledMutationRate(theta_scale);
        thetamax = 0.1;
        if (polydata != nullptr) {
            poissonrandomfield =
                new PoissonRandomField(polydata->GetSampleSizeSet(), *GetCodonStateSpace());
            polyprocess = new PolyProcess(*GetCodonStateSpace(), *polydata, *poissonrandomfield,
                *siteaafitnessarray, *nucmatrix, *theta);
            sitepolysuffstatarray = new PolySuffStatArray(Nsite);
            componentpolysuffstatarray = new PolySuffStatArray(Ncat);
        }

        phyloprocess = new PhyloProcess(
            tree.get(), codondata, branchlength, 0, sitesubmatrixarray, polyprocess);
        phyloprocess->Unfold();

        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatarray = new PathSuffStatArray(Ncat);
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

    //! return current omega value
    double GetOmega() const { return omega; }

    //! return current theta value
    double GetTheta() const { return theta->GetTheta(); }

    //! \brief const access to array of length-pathsuffstats across branches
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! \brief const access to array of suff stats for base mixture
    const DirichletSuffStatArray *GetBaseSuffStatArray() const { return basesuffstatarray; }

    //! \brief const access to array of occupancy suff stats for base mixture
    const OccupancySuffStat *GetBaseOccupancies() const { return baseoccupancy; }

    //-------------------
    // Setting and updating
    // ------------------

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

    //! set omega to new value (multi-gene analyses)
    void SetOmega(double inomega) {
        omega = inomega;
        UpdateCodonMatrices();
    }

    //! set omega hyperparams to new value (multi-gene analyses)
    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape) {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    //! set fraction of positively selected hyperparams to new value (multi-gene analyses)
    void SetDPosOmHyperParameters(double inpi, double inmean, double ininvshape) {
        dposompi = inpi;
        dposomhypermean = inmean;
        dposomhyperinvshape = ininvshape;
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

    //! set nucleotide rates to a new value (multi-gene analyses)
    void SetNucRates(
        const std::vector<double> &innucrelrate, const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

    //! copy nucleotide rates into vectors given as arguments (multi-gene
    //! analyses)
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! set base mixture concentration and center parameters to new value
    //! (multi-gene analyses)
    void SetBaseMixture(const Selector<std::vector<double>> &inbasecenterarray,
        const Selector<double> &inbaseconcentrationarray, const Selector<double> &inbaseweight,
        const Selector<int> &inpermut) {
        basecenterarray->Copy(inbasecenterarray);
        baseconcentrationarray->Copy(inbaseconcentrationarray);
        baseweight->Copy(inbaseweight);
        componentalloc->Permute(inpermut);
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatrices() {
        componentcodonmatrixarray->SetOmega(omega);
        componentcodonmatrixarray->UpdateCodonMatrices();
    }

    //! \brief tell codon matrix k that its parameters have changed and that it
    //! should be updated
    void UpdateCodonMatrix(int k) { (*componentcodonmatrixarray)[k].CorruptMatrix(); }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void Update() {
        UpdateBaseOccupancies();
        UpdateOccupancies();
        UpdateMatrices();
        ResampleSub(1.0);
    }

    void PostPred(std::string name) {
        UpdateBaseOccupancies();
        UpdateOccupancies();
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
        if (omegamode < 2) { total += OmegaLogPrior(); }
        if (polyprocess != nullptr) { total += ThetaLogPrior(); }
        return total;
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    // Branch lengths

    //! log prior over branch lengths
    double BranchLengthsLogPrior() const { return branchlength->GetLogProb(); }

    //! log prior over omega (gamma of mean omegahypermean and inverse shape
    //! omegahyperinvshape)
    double OmegaLogPrior() const {
        double ret = 0;
        if (omegaprior == 0) {
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            ret = alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(omega) -
                  beta * omega;
        } else if (omegaprior == 1) {
            if ((dposompi <= 0) || (dposompi >= 1)) {
                std::cerr << "error in omegalogprior: pi is not 0<pi<1\n";
                exit(1);
            }
            if (omega < 1.0) {
                std::cerr << "error in omegalogprior: omega < 1\n";
                exit(1);
            }
            double dposom = omega - 1.0;
            if (dposom == 0) {
                ret = log(1 - dposompi);
            } else {
                ret = log(dposompi);
                double alpha = 1.0 / dposomhyperinvshape;
                double beta = alpha / dposomhypermean;
                ret += alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(dposom) -
                       beta * dposom;
            }
        } else {
            std::cerr << "error in OmegaLogPrior: unrecognized prior mode\n";
            exit(1);
        }
        return ret;
    }

    //! log prior over theta
    double ThetaLogPrior() const {
        if (theta->GetTheta() > thetamax) {
            return -std::numeric_limits<double>::infinity();
        } else {
            return -log(theta->GetTheta());
        }
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

    //! return log prob only at the tips due to polymorphism of the substitution mapping,
    //! as a function of the current codon substitution process
    double PolySuffStatLogProb() const {
        //! sum over all components to get log prob
        if (polyprocess != nullptr) {
            return componentpolysuffstatarray->GetLogProb(
                *poissonrandomfield, *componentaafitnessarray, *nucmatrix, theta->GetTheta());
        } else {
            return 0;
        }
    }

    //! return log prob only at the tips due to polymorphism of the substitution
    //! mapping, over sites allocated to component k of the mixture
    double ComponentPolySuffStatLogProb(int k) const {
        // sum over all sites allocated to component k
        if (polyprocess != nullptr) {
            return componentpolysuffstatarray->GetVal(k).GetLogProb(
                *poissonrandomfield, componentaafitnessarray->GetVal(k), *nucmatrix, theta->GetTheta());
        } else {
            return 0.0;
        }
    }

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray) +
               PolySuffStatLogProb();
        ;
    }

    //! return log prob of the substitution mappings over sites allocated to
    //! component k of the mixture
    double PathSuffStatLogProb(int k) const {
        return componentpathsuffstatarray->GetVal(k).GetLogProb(
                   componentcodonmatrixarray->GetVal(k)) +
               ComponentPolySuffStatLogProb(k);
        ;
    }

    //! return log prob of first-level mixture components (i.e. all amino-acid
    //! profiles drawn from component k of the base distribution), as a function
    //! of the center and concentration parameters of this component
    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(
                   basecenterarray->GetVal(k), 1.0 / baseconcentrationarray->GetVal(k)) +
               ComponentPolySuffStatLogProb(k);
        ;
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! log prob factor to be recomputed when Theta=4*Ne*u
    double ThetaLogProb() const { return ThetaLogPrior() + PolySuffStatLogProb(); }

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

    //! collect sufficient statistics if substitution mappings across sites
    void CollectSitePathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
        if (polyprocess != nullptr) {
            sitepolysuffstatarray->Clear();
            sitepolysuffstatarray->AddSuffStat(*phyloprocess);
        }
    }

    //! gather site-specific sufficient statistics component-wise
    void CollectComponentPathSuffStat() {
        componentpathsuffstatarray->Clear();
        componentpathsuffstatarray->Add(*sitepathsuffstatarray, *sitealloc);
        if (polyprocess != nullptr) {
            componentpolysuffstatarray->Clear();
            componentpolysuffstatarray->Add(*sitepolysuffstatarray, *sitealloc);
        }
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
        CheckMapping();
    }

    void CheckMapping() const {
        for (int taxon = 0; taxon < Ntaxa; taxon++) {
            for (int site = 0; site < Nsite; site++) {
                int sub_state = phyloprocess->GetPathState(taxon, site);
                int data_state = codondata->GetState(taxon, site);
                std::vector<int> path_state_neighbors =
                    codondata->GetCodonStateSpace()->GetNeighbors(sub_state);
                auto find_data_in_sub_neighbor =
                    find(path_state_neighbors.begin(), path_state_neighbors.end(), data_state);
                bool data_sub_not_neighbors =
                    (find_data_in_sub_neighbor == path_state_neighbors.end());
                if (sub_state != data_state and data_sub_not_neighbors) {
                    std::cerr << "Substitution mapping final state is not even a neighbor of "
                                 "the state given by the alignment"
                              << std::endl;
                }
            }
        }
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            totchrono.Start();
            if (blmode < 2) { MoveBranchLengths(); }

            CollectSitePathSuffStat();
            CollectComponentPathSuffStat();

            if (nucmode < 2) { MoveNucRates(); }

            if (omegamode < 2) { MoveOmega(); }

            if (polyprocess != nullptr) { MoveTheta(); }

            aachrono.Start();
            MoveAAMixture(3);
            aachrono.Stop();

            basechrono.Start();
            if (basemode < 2) { MoveBase(3); }
            basechrono.Stop();

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

    //! MH move on theta
    void MoveTheta() {
        Move::Scaling(theta_scale, 1.0, 10, &AAMutSelDSBDPOmegaModel::ThetaLogProb,
            &AAMutSelDSBDPOmegaModel::NoUpdate, this);
        Move::Scaling(theta_scale, 0.3, 10, &AAMutSelDSBDPOmegaModel::ThetaLogProb,
            &AAMutSelDSBDPOmegaModel::NoUpdate, this);
        assert(theta_scale <= thetamax);
    }

    //! MH move on omega
    void MoveOmega() {
        omegapathsuffstat.Clear();
        omegapathsuffstat.AddSuffStat(*componentcodonmatrixarray, *componentpathsuffstatarray);
        if (omegaprior == 0) {
            GibbsResampleOmega();
        } else {
            MultipleTryMoveOmega(100);
        }
        UpdateCodonMatrices();
    }

    void GibbsResampleOmega() {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(
            alpha + omegapathsuffstat.GetCount(), beta + omegapathsuffstat.GetBeta());
    }

    int MultipleTryMoveOmega(int ntry) {
        /*
        if (omega == 1.0)   {

            // initial prob: (1-pi)*p(D | omega == 1.0)
            // ntry values of dposom
            // final prob: pi* < p(D | omega == 1.0 + dposom) >
        }

        else    {

            // ntry values of dposom, with first value being equal to current omega
        - 1.0
            // initial prob: pi* < p(D | omega == 1.0 + dposom) >
            // final prob: (1-pi)*p(D | omega == 1.0)
        }
        // MH decision rule between omega == 1.0 or omega > 1.0
        // if decision -> omega > 1.0, then randomly sample dposom among ntry values
        and set omega = 1.0 + dposom
        */

        double logp0 = omegapathsuffstat.GetLogProb(1.0);

        double alpha = 1.0 / dposomhyperinvshape;
        double beta = alpha / dposomhypermean;

        std::vector<double> logparray(ntry, 0);
        std::vector<double> dposomarray(ntry, 0);
        double max = 0;
        for (int i = 0; i < ntry; i++) {
            if ((!i) && (omega > 1.0)) {
                dposomarray[i] = omega - 1.0;
            } else {
                dposomarray[i] = Random::Gamma(alpha, beta);
                if (!dposomarray[i]) { dposomarray[i] = 1e-5; }
            }
            logparray[i] = omegapathsuffstat.GetLogProb(1.0 + dposomarray[i]);
            if ((!i) || (max < logparray[i])) { max = logparray[i]; }
        }

        double tot = 0;
        std::vector<double> cumulprob(ntry, 0);
        for (int i = 0; i < ntry; i++) {
            tot += exp(logparray[i] - max);
            cumulprob[i] = tot;
        }
        double logp1 = log(tot / ntry) + max;

        // Multiple-Try Gibbs version
        /*
        double m = (logp0>logp1) ? logp0 : logp1;
        double q0 = (1-dposompi) * exp(logp0-m);
        double q1 = dposompi * exp(logp1-m);
        double p0 = q0 / (q0+q1);

        if (Random::Uniform() < p0) {
            omega = 1.0;
        }
        else    {
            // randomly choose dposom among the ntry values in array
        }
        */

        // Multiple-Try MH version
        int accept = 0;
        int choice = 0;

        if (omega == 1.0) {
            double logratio = log(dposompi) + logp1 - log(1.0 - dposompi) - logp0;
            if (log(Random::Uniform()) < logratio) {
                choice = 1;
                accept = 1;
            } else {
                choice = 0;
            }
        } else {
            double logratio = log(1.0 - dposompi) + logp0 - log(dposompi) - logp1;
            if (log(Random::Uniform()) < logratio) {
                choice = 0;
                accept = 1;
            } else {
                choice = 1;
            }
        }

        if (choice == 1) {
            // randomly choose dposom among the ntry values in array
            double u = tot * Random::Uniform();
            int i = 0;
            while ((i < ntry) && (u > cumulprob[i])) { i++; }
            if (i == ntry) {
                std::cerr << "error in MultipleTryMoveOmega: overflow\n";
                exit(1);
            }
            omega = 1.0 + dposomarray[i];
        } else {
            omega = 1.0;
        }

        return accept;
    }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        Move::Profile(nucrelrate, 0.1, 1, 3, &AAMutSelDSBDPOmegaModel::NucRatesLogProb,
            &AAMutSelDSBDPOmegaModel::UpdateMatrices, this);
        Move::Profile(nucrelrate, 0.03, 3, 3, &AAMutSelDSBDPOmegaModel::NucRatesLogProb,
            &AAMutSelDSBDPOmegaModel::UpdateMatrices, this);
        Move::Profile(nucrelrate, 0.01, 3, 3, &AAMutSelDSBDPOmegaModel::NucRatesLogProb,
            &AAMutSelDSBDPOmegaModel::UpdateMatrices, this);

        Move::Profile(nucstat, 0.1, 1, 3, &AAMutSelDSBDPOmegaModel::NucRatesLogProb,
            &AAMutSelDSBDPOmegaModel::UpdateMatrices, this);
        Move::Profile(nucstat, 0.01, 1, 3, &AAMutSelDSBDPOmegaModel::NucRatesLogProb,
            &AAMutSelDSBDPOmegaModel::UpdateMatrices, this);
    }

    //! MCMC module for the mixture amino-acid fitness profiles
    void MoveAAMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveAAProfiles();
            ResampleEmptyComponents();
            ResampleAlloc();
            LabelSwitchingMove();
            ResampleWeights();
            MoveKappa();
            CollectComponentPathSuffStat();
            UpdateCodonMatrices();
        }
    }

    //! resample empty components of the mixture from prior
    void ResampleEmptyComponents() {
        componentaafitnessarray->PriorResample(*occupancy);
        componentcodonmatrixarray->UpdateCodonMatrices(*occupancy);
    }

    //! MH move on amino-acid fitness profiles (occupied components only)
    void MoveAAProfiles() {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double CompMoveAAProfiles(int nrep) {
        accb1 += MoveAA(1.0, 1, nrep);
        // accb2 += MoveAA(1.0,3,nrep);
        // accb3 += MoveAA(0.3,3,nrep);
        accb4 += MoveAA(0.1, 3, nrep);
        totb1++;
        totb2++;
        totb3++;
        totb4++;
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MulMoveAAProfiles(int nrep) {
        acca1 += MoveAAGamma(3.0, nrep);
        acca2 += MoveAAGamma(1.0, nrep);
        // acca3 += MoveAAGamma(0.3,nrep);
        // acca4 += MoveAAGamma(0.1,nrep);
        tota1++;
        tota2++;
        tota3++;
        tota4++;
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double MoveAA(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Naa];
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                std::vector<double> &aa = (*componentaafitnessarray)[i];
                for (int rep = 0; rep < nrep; rep++) {
                    for (int l = 0; l < Naa; l++) { bk[l] = aa[l]; }
                    double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                    deltalogprob += loghastings;
                    UpdateCodonMatrix(i);
                    deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        for (int l = 0; l < Naa; l++) { aa[l] = bk[l]; }
                        UpdateCodonMatrix(i);
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
            if (occupancy->GetVal(i)) {
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

                    UpdateCodonMatrix(i);

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
                        UpdateCodonMatrix(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! Gibbs resample mixture allocations
    void ResampleAlloc() {
        std::vector<double> postprob(Ncat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetAllocPostProb(i, postprob);
            sitealloc->GibbsResample(i, postprob);
        }
        UpdateOccupancies();
    }

    //! update mixture occupancy suff stats (for resampling mixture weights)
    void UpdateOccupancies() {
        occupancy->Clear();
        occupancy->AddSuffStat(*sitealloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetAllocPostProb(int site, std::vector<double> &postprob) {
        double max = 0;
        const std::vector<double> &w = weight->GetArray();
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        for (int i = 0; i < Ncat; i++) {
            double tmp = suffstat.GetLogProb(componentcodonmatrixarray->GetVal(i));
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
    void LabelSwitchingMove() {
        Permutation permut(Ncat);
        weight->LabelSwitchingMove(5, *occupancy, permut);
        sitealloc->Permute(permut);
        componentaafitnessarray->Permute(permut);
    }

    //! Gibbs resample mixture weights (based on occupancy suff stats)
    void ResampleWeights() { weight->GibbsResample(*occupancy); }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa() {
        Move::Scaling(kappa, 1.0, 10, &AAMutSelDSBDPOmegaModel::StickBreakingHyperLogProb,
            &AAMutSelDSBDPOmegaModel::NoUpdate, this);
        Move::Scaling(kappa, 0.3, 10, &AAMutSelDSBDPOmegaModel::StickBreakingHyperLogProb,
            &AAMutSelDSBDPOmegaModel::NoUpdate, this);
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
            if (baseoccupancy->GetVal(k)) {
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
            if (baseoccupancy->GetVal(k)) {
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
        basecenterarray->PriorResample(*baseoccupancy);
        baseconcentrationarray->PriorResample(*baseoccupancy);
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

    //! update base occupancy suff stats (for moving base weights)
    void UpdateBaseOccupancies() {
        baseoccupancy->Clear();
        baseoccupancy->AddSuffStat(*componentalloc);
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
        baseweight->LabelSwitchingMove(5, *baseoccupancy, permut);
        componentalloc->Permute(permut);
        basecenterarray->Permute(permut);
        baseconcentrationarray->Permute(permut);
        basesuffstatarray->Permute(permut);
    }

    //! Gibbs resample base mixture weights (based on occupancy suff stats)
    void ResampleBaseWeights() { baseweight->GibbsResample(*baseoccupancy); }

    //! MH move on basekappa, concentration parameter of the base mixture
    void MoveBaseKappa() {
        Move::Scaling(basekappa, 1.0, 10, &AAMutSelDSBDPOmegaModel::BaseStickBreakingHyperLogProb,
            &AAMutSelDSBDPOmegaModel::NoUpdate, this);
        Move::Scaling(basekappa, 0.3, 10, &AAMutSelDSBDPOmegaModel::BaseStickBreakingHyperLogProb,
            &AAMutSelDSBDPOmegaModel::NoUpdate, this);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    //! return number of occupied components in first-level mixture (mixture of
    //! amino-acid fitness profiles)
    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) { n++; }
        }
        return n;
    }

    //! return number of occupied components in base distribution
    int GetBaseNcluster() const {
        int n = 0;
        for (int i = 0; i < baseNcat; i++) {
            if (baseoccupancy->GetVal(i)) { n++; }
        }
        return n;
    }

    //! return mean entropy of amino-acd fitness profiles
    double GetMeanAAEntropy() const { return componentaafitnessarray->GetMeanEntropy(); }

    //! return mean of concentration parameters of base distribution
    double GetMeanComponentAAConcentration() const {
        double tot = 0;
        double totw = 0;
        for (int i = basemin; i < baseNcat; i++) {
            tot += baseoccupancy->GetVal(i) * baseconcentrationarray->GetVal(i);
            totw += baseoccupancy->GetVal(i);
        }
        return tot / totw;
    }

    //! return mean entropy of centers of base distribution
    double GetMeanComponentAAEntropy() const {
        double tot = 0;
        for (int i = 0; i < baseNcat; i++) {
            tot += baseoccupancy->GetVal(i) * Random::GetEntropy(basecenterarray->GetVal(i));
        }
        return tot / Ncat;
    }

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const { return Random::GetEntropy(nucrelrate); }

    double GetPredictedDNDS() const {
        double mean = 0;
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                mean += occupancy->GetVal(i) * (*componentcodonmatrixarray)[i].GetPredictedDNDS();
            }
        }
        mean /= Nsite;
        return mean;
    }

    const std::vector<double> &GetProfile(int i) const { return siteaafitnessarray->GetVal(i); }

    void Monitor(std::ostream &os) const {
        os << totchrono.GetTime() << '\t' << aachrono.GetTime() << '\t' << basechrono.GetTime()
           << '\n';
        os << "prop time in aa moves  : " << aachrono.GetTime() / totchrono.GetTime() << '\n';
        os << "prop time in base moves: " << basechrono.GetTime() / totchrono.GetTime() << '\n';
    }

    void ToStream(std::ostream &os) const {
        os << "AAMutSelDSBDPOmega" << '\t';
        os << datafile << '\t' << treefile << '\t';
        os << omegamode << '\t' << omegaprior << '\t';
        os << dposompi << '\t' << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
        os << Ncat << '\t' << baseNcat << '\t';
        os << polymorphism_aware << '\t';
        tracer->write_line(os);
    }

    AAMutSelDSBDPOmegaModel(std::istream &is) {
        std::string model_name;
        is >> model_name;
        if (model_name != "AAMutSelDSBDPOmega") {
            std::cerr << "Expected AAMutSelDSBDPOmega for model name, got " << model_name << "\n";
            exit(1);
        }
        is >> datafile >> treefile;
        is >> omegamode >> omegaprior;
        is >> dposompi >> dposomhypermean >> dposomhyperinvshape;
        is >> Ncat >> baseNcat;
        is >> polymorphism_aware;

        init();
        tracer->read_line(is);
        Update();
    };
};
