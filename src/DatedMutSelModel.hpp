#pragma once

#include <assert.h>
#include "AAMutSelNeCodonMatrixBidimArray.hpp"
#include "BranchComponentMatrixSelector.hpp"
#include "BranchProduct.hpp"
#include "Chrono.hpp"
#include "Chronogram.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "Move.hpp"
#include "MultinomialAllocationVector.hpp"
#include "NodeMultivariateProcess.hpp"
#include "Permutation.hpp"
#include "PhyloProcess.hpp"
#include "PolyProcess.hpp"
#include "PolySuffStat.hpp"
#include "ScaledMutationRateCompound.hpp"
#include "ScatterSuffStat.hpp"
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
 * - A branch Ne
 * - A branch mutation rate (correlated to Ne)
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
 */

class DatedMutSelModel : public ChainComponent {
    std::string datafile, treefile;

    bool condition_aware;
    bool polymorphism_aware;

    std::unique_ptr<const Tree> tree;

    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;
    PolyData *polydata{nullptr};

    int Nsite;
    int Ntaxa;
    int Nbranch;
    // number of conditions (each with different Ne)

    // Node ages
    NodeAges *nodeages;
    // Chronogram (diff between node ages)
    Chronogram *chronogram;

    int dimension;
    int invert_whishart_df;
    double invert_whishart_kappa;
    // Covariance matrix
    EMatrix precision_matrix;
    EVector root_multivariate;
    NodeMultivariateProcess *node_multivariate;

    // Branch rates (brownian process)
    NodeProcess *noderates;
    BranchProcess *branchrates;

    // Branch lengths (product of branch rates and chronogram)
    BranchwiseProduct *branchlength;
    BranchwiseProduct *branchscaledpopsize;

    // Branch Population size (brownian process)
    NodeProcess *nodepopsize;
    BranchProcess *branchpopsize;

    // nucleotide rates hyperparameters
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucstat;
    std::vector<double> nucrelrate;
    GTRSubMatrix *nucmatrix;

    // base distribution G0 is itself a stick-breaking mixture of Dirichlet
    // distributions
    int baseNcat;
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

    // which site is under which component
    MultinomialAllocationVector *sitealloc;

    // Bi-dimensional array of codon matrices (one for each distinct branch condition, and one for
    // each aa fitness profile)
    MutSelNeCodonMatrixBidimArray *branchcomponentcodonmatrixarray;
    AAMutSelNeCodonSubMatrixArray *rootcomponentcodonmatrixarray;

    // this one is used by PhyloProcess: has to be a BranchComponentMatrixSelector<SubMatrix>
    BranchComponentMatrixSelector<SubMatrix> *branchsitecodonmatrixarray;

    // this one is also used by PhyloProcess: has to be a RootComponentMatrixSelector<SubMatrix>
    RootComponentMatrixSelector<SubMatrix> *rootsitecodonmatrixarray;

    // this one is used by PolyProcess: has to be a Selector<vector<double>>
    MixtureSelector<std::vector<double>> *siteaafitnessarray;

    PhyloProcess *phyloprocess;

    // global theta (4*Ne*u) used for polymorphism
    double theta_scale;
    double thetamax;
    CompoundScaledMutationRate *theta;

    PolyProcess *polyprocess{nullptr};
    PoissonRandomField *poissonrandomfield{nullptr};

    PathSuffStatBidimArray *sitepathsuffstatbidimarray;
    PathSuffStatBidimArray *componentpathsuffstatbidimarray;

    PolySuffStatBidimArray *sitepolysuffstatbidimarray{nullptr};
    PolySuffStatBidimArray *componentpolysuffstatbidimarray{nullptr};

    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    ScatterSuffStat *scattersuffstat;

  public:
    friend std::ostream &operator<<(std::ostream &os, DatedMutSelModel &m);

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
    DatedMutSelModel(std::string indatafile, std::string intreefile, int inNcat, int inbaseNcat,
        bool incondition_aware, bool inpolymorphism_aware)
        : datafile(indatafile),
          treefile(intreefile),
          condition_aware(incondition_aware),
          polymorphism_aware(inpolymorphism_aware),
          baseNcat(inbaseNcat),
          Ncat(inNcat) {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        if (polymorphism_aware) { polydata = new PolyData(codondata, datafile); }

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        if (Ncat <= 0) { Ncat = Nsite; }
        if (Ncat > Nsite) { Ncat = Nsite; }
        if (Ncat > 100) { Ncat = 100; }

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        std::cerr << "ncat : " << Ncat << '\n';
        std::cerr << "basencat : " << baseNcat << '\n';

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
        invert_whishart_df = dimension + 1;
        invert_whishart_kappa = 1.0;
        root_multivariate = EVector::Constant(dimension, -1.0);
        root_multivariate(1) = 0.0;
        precision_matrix = EMatrix::Identity(dimension, dimension) * 10;

        node_multivariate =
            new NodeMultivariateProcess(*chronogram, precision_matrix, root_multivariate);

        // Branch rates (brownian process)
        noderates = new NodeProcess(*node_multivariate, 0);
        branchrates = new BranchProcess(*noderates);

        // Branch lengths (product of branch rates and chronogram)
        branchlength = new BranchwiseProduct(*chronogram, *branchrates);

        // Branch omega (brownian process)
        nodepopsize = new NodeProcess(*node_multivariate, 1);
        branchpopsize = new BranchProcess(*nodepopsize);
        branchscaledpopsize = new BranchwiseProduct(*chronogram, *branchpopsize);

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

        // codon matrices per branch and per site
        branchcomponentcodonmatrixarray = new MutSelNeCodonMatrixBidimArray(
            GetCodonStateSpace(), nucmatrix, componentaafitnessarray, &branchpopsize->GetArray());

        // sub matrices per branch and per site
        branchsitecodonmatrixarray = new BranchComponentMatrixSelector<SubMatrix>(
            branchcomponentcodonmatrixarray, sitealloc, *tree);

        rootcomponentcodonmatrixarray = new AAMutSelNeCodonSubMatrixArray(GetCodonStateSpace(),
            nucmatrix, componentaafitnessarray, nodepopsize->GetExpVal(tree->root()));

        // sub matrices for root, across sites
        rootsitecodonmatrixarray =
            new RootComponentMatrixSelector<SubMatrix>(rootcomponentcodonmatrixarray, sitealloc);

        // selector, specifying which aa fitness array should be used for each site
        siteaafitnessarray =
            new MixtureSelector<std::vector<double>>(componentaafitnessarray, sitealloc);

        // global theta (4*Ne*u = 1e-5 by default, and maximum value 0.1)
        theta_scale = 1e-5;
        theta = new CompoundScaledMutationRate(Ntaxa, theta_scale, noderates, nodepopsize);
        thetamax = 0.1;
        if (polydata != nullptr) {
            poissonrandomfield =
                new PoissonRandomField(polydata->GetSampleSizeSet(), *GetCodonStateSpace());
            polyprocess = new PolyProcess(*GetCodonStateSpace(), *polydata, *poissonrandomfield,
                *siteaafitnessarray, *nucmatrix, *theta);
            sitepolysuffstatbidimarray = new PolySuffStatBidimArray(Nsite, Ntaxa);
            componentpolysuffstatbidimarray = new PolySuffStatBidimArray(Ncat, Ntaxa);
        }

        phyloprocess = new PhyloProcess(tree.get(), codondata, branchlength, nullptr,
            branchsitecodonmatrixarray, rootsitecodonmatrixarray, polyprocess);
        phyloprocess->Unfold();

        sitepathsuffstatbidimarray = new PathSuffStatBidimArray(Nbranch, Nsite);
        componentpathsuffstatbidimarray = new PathSuffStatBidimArray(Nbranch, Ncat);
        scattersuffstat = new ScatterSuffStat(*tree);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
    }

    virtual ~DatedMutSelModel() = default;

    void move(int it) override { Move(); }

    template <class C>
    void declare_model(C &t) {
        t.add("nodeages", *nodeages);
        t.add("node_multivariate", *node_multivariate);
        t.add("covmatrix", precision_matrix);
        t.add("root_multivariate", root_multivariate);
        t.add("invert_whishart_kappa", invert_whishart_kappa);
        t.add("invert_whishart_df", invert_whishart_df);
        t.add("nucrelrate", nucrelrate);
        t.add("nucstat", nucstat);
        t.add("basekappa", basekappa);
        t.add("baseweight", *baseweight);
        t.add("componentalloc", *componentalloc);
        t.add("*basecenterarray", *basecenterarray);
        t.add("*baseconcentrationarray", *baseconcentrationarray);
        t.add("kappa", kappa);
        t.add("weight", *weight);
        t.add("componentaafitnessarray", *componentaafitnessarray);
        t.add("sitealloc", *sitealloc);
        if (polyprocess != nullptr) { t.add("theta_scale", theta_scale); }
    }

    template <class C>
    void declare_stats(C &t) {
        t.add("logprior", [this]() { return GetLogPrior(); });
        t.add("lnL", [this]() { return GetLogLikelihood(); });
        // 3x: per coding site (and not per nucleotide site)
        t.add("blengthmean", [&]() { return branchlength->GetMean(); });
        t.add("blengthvar", [&]() { return branchlength->GetVar(); });
        t.add("predicted_dnds", [this]() { return GetPredictedDNDS(); });
        if (polyprocess != nullptr) { t.add("theta_scale", theta_scale); }
        t.add("ncluster", [this]() { return GetNcluster(); });
        t.add("kappa", kappa);
        if (baseNcat > 1) {
            t.add("basencluster", [this]() { return GetBaseNcluster(); });
            t.add("basekappa", basekappa);
        }
        for (int i = 0; i < dimension; i++) {
            t.add("root_" + std::to_string(i), root_multivariate.coeffRef(i));
            for (int j = 0; j <= i; j++) {
                t.add("cov_" + std::to_string(i) + "_" + std::to_string(j),
                    precision_matrix.coeffRef(i, j));
            }
        }
        t.add("bomegamean", [&]() { return branchscaledpopsize->GetMean(); });
        t.add("bomegavar", [&]() { return branchscaledpopsize->GetVar(); });
        t.add("aaent", [this]() { return GetMeanAAEntropy(); });
        t.add("meanaaconc", [this]() { return GetMeanComponentAAConcentration(); });
        t.add("aacenterent", [this]() { return GetMeanComponentAAEntropy(); });
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

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //-------------------
    // Setting and updating
    // ------------------

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatrices() {
        branchcomponentcodonmatrixarray->SetNe(branchpopsize->GetArray());
        branchcomponentcodonmatrixarray->CorruptCodonMatrices();
        rootcomponentcodonmatrixarray->SetNe(nodepopsize->GetExpVal(tree->root()));
        rootcomponentcodonmatrixarray->UpdateCodonMatrices();
    }

    //! \brief tell codon matrices for site i and across conditions that their parameters have
    //! changed and that they should be updated
    void UpdateCodonMatrix(int i) {
        branchcomponentcodonmatrixarray->CorruptColCodonMatrices(i);
        rootcomponentcodonmatrixarray->UpdateCodonMatrices(i);
    }

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

    //! \brief Update the BranchProcess (rates and omega) with the underlying NodeProcess, Update
    //! the Chronogram with the underlying NodeAges. And finally update the branch lengths with the
    //! Chronogram and the BranchProcess (rates).
    //!
    //! Used when the model is restarted or for the posterior predictif.
    void UpdateBranches() {
        chronogram->Update();
        branchpopsize->Update();
        branchscaledpopsize->Update();
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
    void UpdateLocalBranchPopSize(Tree::NodeIndex node) {
        branchpopsize->UpdateLocal(node);
        if (!tree->is_root(node)) {
            Tree::BranchIndex branch = tree->branch_index(node);
            branchcomponentcodonmatrixarray->SetRowNe(branch, branchpopsize->GetVal(branch));
            branchcomponentcodonmatrixarray->CorruptRowCodonMatrices(branch);
        } else {
            rootcomponentcodonmatrixarray->SetNe(nodepopsize->GetExpVal(tree->root()));
            rootcomponentcodonmatrixarray->UpdateCodonMatrices();
        }
    }

    void Update() {
        UpdateBranches();
        UpdateBaseOccupancies();
        UpdateOccupancies();
        UpdateMatrices();
        ResampleSub(1.0);
    }

    void PostPred(std::string name) {
        UpdateBranches();
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
        total += NodeMultivariateLogPrior();
        total += RootMultivariateLogPrior();
        total += NucRatesLogPrior();

        if (baseNcat > 1) {
            total += BaseStickBreakingHyperLogPrior();
            total += BaseStickBreakingLogPrior();
        }
        total += BaseLogPrior();
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AALogPrior();

        if (polyprocess != nullptr) { total += ThetaLogPrior(); }
        return total;
    }

    //! return log likelihood
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
    double RootMultivariateHyperLogPrior() const { return -root_multivariate.sum(); }


    //! log prior over theta
    double ThetaLogPrior() const {
        if (theta_scale > thetamax) {
            return -std::numeric_limits<double>::infinity();
        } else {
            return -log(theta_scale);
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
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics if substitution mappings across sites
    void CollectSitePathSuffStat() {
        sitepathsuffstatbidimarray->Clear();
        sitepathsuffstatbidimarray->AddSuffStat(*phyloprocess);
        if (polyprocess != nullptr) {
            sitepolysuffstatbidimarray->Clear();
            sitepolysuffstatbidimarray->AddSuffStat(*phyloprocess);
        }
    }

    //! gather site-specific sufficient statistics component-wise
    void CollectComponentPathSuffStat() {
        componentpathsuffstatbidimarray->Clear();
        componentpathsuffstatbidimarray->Add(*sitepathsuffstatbidimarray, *sitealloc);
        if (polyprocess != nullptr) {
            componentpolysuffstatbidimarray->Clear();
            componentpolysuffstatbidimarray->Add(*sitepolysuffstatbidimarray, *sitealloc);
        }
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    // Scatter (brownian process)

    void CollectScatterSuffStat() {
        scattersuffstat->Clear();
        scattersuffstat->AddSuffStat(*node_multivariate);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! return log prob only at the tips due to polymorphism of the substitution mapping,
    //! as a function of the current codon substitution process
    double PolySuffStatLogProb() const {
        //! sum over all components to get log prob
        if (polyprocess != nullptr) {
            return componentpolysuffstatbidimarray->GetLogProb(
                *poissonrandomfield, *componentaafitnessarray, *nucmatrix, *theta);
        } else {
            return 0;
        }
    }

    //! return log prob only at the tips due to polymorphism of the substitution
    //! mapping, over sites allocated to component k of the mixture
    double ComponentPolySuffStatLogProb(int k) const {
        // sum over all sites allocated to component k
        if (polyprocess != nullptr) {
            double tot = 0;
            for (int taxon = 0; taxon < Ntaxa; taxon++) {
                tot += componentpolysuffstatbidimarray->GetVal(k, taxon).GetLogProb(
                    *poissonrandomfield, componentaafitnessarray->GetVal(k), *nucmatrix,
                    theta->GetTheta(taxon));
            }
            return tot;
        } else {
            return 0.0;
        }
    }

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return componentpathsuffstatbidimarray->GetLogProb(*branchcomponentcodonmatrixarray) +
               PolySuffStatLogProb();
    }

    //! return log prob of the substitution mappings over sites allocated to
    //! component k of the mixture
    double ComponentPathSuffStatLogProb(int k) const {
        double loglk = 0.0;
        for (int cond{0}; cond < Nbranch; cond++) {
            loglk += componentpathsuffstatbidimarray->GetVal(cond, k).GetLogProb(
                branchcomponentcodonmatrixarray->GetVal(cond, k));
        }
        loglk += ComponentPolySuffStatLogProb(k);
        return loglk;
    }

    //! return log prob of first-level mixture components (i.e. all amino-acid
    //! profiles drawn from component k of the base distribution), as a function
    //! of the center and concentration parameters of this component
    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(
                   basecenterarray->GetVal(k), 1.0 / baseconcentrationarray->GetVal(k)) +
               ComponentPolySuffStatLogProb(k);
    }

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

    // PopSize
    //! \brief log prob to be recomputed when moving omega (brownian process) around of focal node
    double LocalNodePopSizeLogProb(Tree::NodeIndex node) const {
        return LocalNodeMultivariateLogPrior(node) + LocalNodePopSizeSuffStatLogProb(node);
    }

    //! \brief log prob factor (without prior) to be recomputed when moving omega (brownian process)
    //! around of focal node
    double LocalNodePopSizeSuffStatLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        // for all children
        for (auto const &child : tree->children(node)) {
            tot += NodePopSizeSuffStatLogProb(tree->branch_index(child));
        }
        if (!tree->is_root(node)) {
            // for the branch attached to the node
            tot += NodePopSizeSuffStatLogProb(tree->branch_index(node));
        }
        assert(tot != 0);
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! omega of a given branch
    double NodePopSizeSuffStatLogProb(Tree::BranchIndex branch) const {
        return componentpathsuffstatbidimarray->GetRowLogProb(
            branch, *branchcomponentcodonmatrixarray);
    }

    //! Root log prob
    double RootLogProb() const {
        return RootMultivariateHyperLogPrior() + RootMultivariateLogPrior();
    }


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
            CollectLengthSuffStat();
            MoveNodeAges(1.0, 3);
            MoveNodeAges(0.1, 3);

            MoveNodeRates(1.0, 3);
            MoveNodeRates(0.1, 3);

            CollectSitePathSuffStat();
            CollectComponentPathSuffStat();
            MoveNucRates();
            MoveNodePopSize(1.0, 3);
            MoveNodePopSize(0.1, 3);

            CollectScatterSuffStat();
            SamplePrecisionMatrix();
            MoveRootMultivariate(1.0, 3);
            MoveRootMultivariate(0.1, 3);
            UpdateMatrices();

            if (polyprocess != nullptr) { MoveTheta(); }

            MoveAAMixture(3);
            MoveBase(3);
        }
    }


    void SamplePrecisionMatrix() {
        scattersuffstat->SamplePrecisionMatrix(
            precision_matrix, invert_whishart_df, invert_whishart_kappa);
    };

    void MoveRootMultivariate(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (int dim = 0; dim < dimension; dim++) {
                double logratio = -RootLogProb();

                double m = tuning * (Random::Uniform() - 0.5);

                root_multivariate(dim) += m;

                logratio += RootLogProb();

                bool accept = (log(Random::Uniform()) < logratio);
                if (!accept) { root_multivariate(dim) -= m; }
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
        double logratio = -LocalNodeAgeLogProb(node);

        double bk = nodeages->GetVal(node);
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
            noderates->SlidingMove(node, -m);
            UpdateLocalBranchRates(node);
        }
    }

    //! MH moves on branch omega (brownian process)
    void MoveNodePopSize(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(tree->nb_nodes()); node++) {
                MoveNodePopSize(node, tuning);
            }
        }
    }

    //! MH moves on branch omega (brownian process) for a focal node
    void MoveNodePopSize(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodePopSizeLogProb(node);

        double m = tuning * nodepopsize->GetSigma() * (Random::Uniform() - 0.5);
        nodepopsize->SlidingMove(node, m);
        UpdateLocalBranchPopSize(node);

        logratio += LocalNodePopSizeLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            nodepopsize->SlidingMove(node, -m);
            UpdateLocalBranchPopSize(node);
        }
    }

    //! MH move on base mixture
    void MoveBase(int nrep) {
        if (baseNcat > 1) { ResampleBaseAlloc(); }
        MoveBaseMixture(nrep);
    }


    //! MH move on theta
    void MoveTheta() {
        Move::Scaling(theta_scale, 1.0, 10, &DatedMutSelModel::ThetaLogProb,
            &DatedMutSelModel::NoUpdate, this);
        Move::Scaling(theta_scale, 0.3, 10, &DatedMutSelModel::ThetaLogProb,
            &DatedMutSelModel::NoUpdate, this);
        assert(theta_scale <= thetamax);
    }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        Move::Profile(nucrelrate, 0.1, 1, 3, &DatedMutSelModel::NucRatesLogProb,
            &DatedMutSelModel::UpdateMatrices, this);
        Move::Profile(nucrelrate, 0.03, 3, 3, &DatedMutSelModel::NucRatesLogProb,
            &DatedMutSelModel::UpdateMatrices, this);
        Move::Profile(nucrelrate, 0.01, 3, 3, &DatedMutSelModel::NucRatesLogProb,
            &DatedMutSelModel::UpdateMatrices, this);

        Move::Profile(nucstat, 0.1, 1, 3, &DatedMutSelModel::NucRatesLogProb,
            &DatedMutSelModel::UpdateMatrices, this);
        Move::Profile(nucstat, 0.01, 1, 3, &DatedMutSelModel::NucRatesLogProb,
            &DatedMutSelModel::UpdateMatrices, this);
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
        branchcomponentcodonmatrixarray->UpdateCodonMatrices(*occupancy);
    }

    //! MH move on amino-acid fitness profiles (occupied components only)
    void MoveAAProfiles() {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    void CompMoveAAProfiles(int nrep) {
        MoveAA(1.0, 1, nrep);
        MoveAA(0.1, 3, nrep);
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    void MulMoveAAProfiles(int nrep) {
        MoveAAGamma(3.0, nrep);
        MoveAAGamma(1.0, nrep);
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
                    double deltalogprob = -AALogPrior(i) - ComponentPathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                    deltalogprob += loghastings;
                    UpdateCodonMatrix(i);
                    deltalogprob += AALogPrior(i) + ComponentPathSuffStatLogProb(i);
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
                        -GammaAALogPrior(x, aacenter, aaconc) - ComponentPathSuffStatLogProb(i);

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

                    deltalogprob +=
                        GammaAALogPrior(x, aacenter, aaconc) + ComponentPathSuffStatLogProb(i);

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

        // !!!!! Here the branch should not matter, so we use the root.
        const PathSuffStat &suffstat = sitepathsuffstatbidimarray->GetVal(0, site);
        for (int i = 0; i < Ncat; i++) {
            // !!!!! Here the condition should not matter, so we use site 0.
            double tmp = suffstat.GetLogProb(branchcomponentcodonmatrixarray->GetVal(0, 0));
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
        Move::Scaling(kappa, 1.0, 10, &DatedMutSelModel::StickBreakingHyperLogProb,
            &DatedMutSelModel::NoUpdate, this);
        Move::Scaling(kappa, 0.3, 10, &DatedMutSelModel::StickBreakingHyperLogProb,
            &DatedMutSelModel::NoUpdate, this);
    }

    //! MCMC module for the base mixture
    void MoveBaseMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseComponents(10);
            ResampleBaseEmptyComponents();
            if (baseNcat > 1) {
                BaseLabelSwitchingMove();
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
        for (int k = 0; k < baseNcat; k++) {
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
        Move::Scaling(basekappa, 1.0, 10, &DatedMutSelModel::BaseStickBreakingHyperLogProb,
            &DatedMutSelModel::NoUpdate, this);
        Move::Scaling(basekappa, 0.3, 10, &DatedMutSelModel::BaseStickBreakingHyperLogProb,
            &DatedMutSelModel::NoUpdate, this);
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
        for (int i = 0; i < baseNcat; i++) {
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

    double GetPredictedDNDS() const {
        double mean = 0;
        for (int k = 0; k < Ncat; k++) {
            if (occupancy->GetVal(k)) {
                for (int cond{0}; cond < Nbranch; cond++) {
                    mean += occupancy->GetVal(k) *
                            branchcomponentcodonmatrixarray->GetVal(cond, k).GetPredictedDNDS();
                }
            }
        }
        mean /= (Nsite * Nbranch);
        assert(mean <= 1);
        return mean;
    }

    const std::vector<double> &GetProfile(int i) const { return siteaafitnessarray->GetVal(i); }

    void ToStream(std::ostream &os) { os << *this; }
};

std::istream &operator>>(std::istream &is, std::unique_ptr<DatedMutSelModel> &m) {
    std::string model_name, datafile, treefile;
    int Ncat, baseNcat;
    bool condition_aware, polymorphism_aware;

    is >> model_name;
    if (model_name != "DatedMutSelModel") {
        std::cerr << "Expected DatedMutSelModel for model name, got " << model_name << "\n";
        exit(1);
    }

    is >> datafile >> treefile;
    is >> Ncat >> baseNcat;
    is >> condition_aware >> polymorphism_aware;
    m.reset(new DatedMutSelModel(
        datafile, treefile, Ncat, baseNcat, condition_aware, polymorphism_aware));
    Tracer tracer{*m, &DatedMutSelModel::declare_model};
    tracer.read_line(is);
    m->Update();
    return is;
}

std::ostream &operator<<(std::ostream &os, DatedMutSelModel &m) {
    Tracer tracer{m, &DatedMutSelModel::declare_model};
    os << "DatedMutSelModel" << '\t';
    os << m.datafile << '\t';
    os << m.treefile << '\t';
    os << m.Ncat << '\t';
    os << m.baseNcat << '\t';
    os << m.condition_aware << '\t';
    os << m.polymorphism_aware << '\t';
    tracer.write_line(os);
    return os;
}
