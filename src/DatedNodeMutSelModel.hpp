#pragma once

#include <cassert>
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
#include "MultivariateProcess.hpp"
#include "Permutation.hpp"
#include "PhyloProcess.hpp"
#include "PolyProcess.hpp"
#include "PolySuffStat.hpp"
#include "ScaledMutationRateCompound.hpp"
#include "ScatterSuffStat.hpp"
#include "StickBreakingProcess.hpp"
#include "TaxonTraits.hpp"
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

double distance(std::vector<double> const &v1, std::vector<double> const &v2) {
    double tot = 0;
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < v1.size(); i++) { tot += abs(v1[i] - v2[i]); }
    return tot;
}

std::tuple<std::vector<std::vector<double>>, std::vector<size_t>> open_preferences(
    std::string const &file_name) {
    std::vector<std::vector<double>> fitness_profiles{};
    std::vector<size_t> alloc;

    std::ifstream input_stream(file_name);
    if (!input_stream) {
        std::cerr << "Preferences file " << file_name << " doesn't exist" << std::endl;
    }

    std::string line;

    // skip the header of the file
    getline(input_stream, line);
    char sep{' '};
    long nbr_col = 0;
    for (char sep_test : std::vector<char>({' ', ',', '\t'})) {
        long n = std::count(line.begin(), line.end(), sep_test);
        if (n > nbr_col) {
            sep = sep_test;
            nbr_col = n + 1;
        }
    }
    nbr_col -= 20;

    while (getline(input_stream, line)) {
        std::vector<double> fitness_profil(20, 0.0);
        std::string word;
        std::istringstream line_stream(line);
        unsigned counter{0};

        while (getline(line_stream, word, sep)) {
            if (counter >= nbr_col) { fitness_profil[counter - nbr_col] = stod(word); }
            counter++;
        }

        bool push = true;
        assert(std::abs(std::accumulate(fitness_profil.begin(), fitness_profil.end(), 0.0) - 1.0) <
               1e-5);
        for (size_t i = 0; i < fitness_profiles.size(); i++) {
            if (distance(fitness_profiles[i], fitness_profil) < 1e-5) {
                push = false;
                alloc.push_back(i);
                break;
            }
        }
        if (push) {
            alloc.push_back(fitness_profiles.size());
            fitness_profiles.push_back(fitness_profil);
        }
        assert(alloc[alloc.size() - 1] < fitness_profiles.size());
    }
    return std::make_tuple(fitness_profiles, alloc);
}

class DatedNodeMutSelModel : public ChainComponent {
    std::string datafile, treefile, traitsfile{"Null"}, profilesfile{"Null"}, fossilsfile{"Null"};

    bool condition_aware;
    bool polymorphism_aware;
    unsigned precision;
    bool arithmetic;
    bool clamp_profiles{false}, move_root_pop_size, clamp_pop_sizes, clamp_nuc_matrix,
        clamp_corr_matrix;
    std::unique_ptr<const Tree> tree;

    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;
    TaxonTraits *taxon_traits{nullptr};
    PolyData *polydata{nullptr};

    int Nsite;
    int Ntaxa;
    int Nbranch;
    // number of conditions (each with different Ne)

    // Node ages
    NodeAges *nodeages;
    // Chronogram (diff between node ages)
    Chronogram *chronogram;

    // Precision matrix and prior
    int dimensions;
    int prior_cov_df;
    bool uniq_kappa;
    PriorCovariance *prior_matrix;
    PrecisionMatrix *precision_matrix;
    NodeMultivariateProcess *node_multivariate;

    // Branch Population size (brownian process)
    NodeProcess *nodepopsize;
    BranchProcess *branchpopsize;
    double root_popsize{1.0};

    // Branch mutation rates (nbr of mutations per generation)
    NodeProcess *nodemutrates;
    BranchProcess *branchmutrates;

    // Branch generation rate (nbr of generations per time)
    NodeProcess *nodegentimes;
    BranchProcess *branchgentimes;

    // Branch lengths (product of branch rates and chronogram)
    BranchwiseProduct *branchlength;

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
    NodeProcessScaledMutationRate *theta;

    PolyProcess *polyprocess{nullptr};
    PoissonRandomField *poissonrandomfield{nullptr};

    PathSuffStatBidimArray *branchcomponentpathsuffstatbidimarray;
    PathSuffStatBidimArray *branchsitepathsuffstatbidimarray;

    PathSuffStatArray *rootcomponentpathsuffstatarray;
    PathSuffStatArray *rootsitepathsuffstatarray;

    PolySuffStatBidimArray *taxoncomponentpolysuffstatbidimarray{nullptr};
    PolySuffStatBidimArray *taxonsitepolysuffstatbidimarray{nullptr};

    PoissonSuffStatBranchArray *branchlengthpathsuffstatarray;
    ScatterSuffStat *scattersuffstat;

    SimpleBranchArray<double> *branchdnds;

    std::map<std::string, double> moves_tot;
    std::map<std::string, double> moves_accepted;
    std::map<std::string, Chrono> moves_chrono;
    std::map<std::string, double> moves_chrono_time;
    double chronoresamblesub;
    double chronomove;

  public:
    friend std::ostream &operator<<(std::ostream &os, DatedNodeMutSelModel &m);

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
    DatedNodeMutSelModel(std::string indatafile, std::string intreefile, std::string intraitsfile,
        std::string inprofiles, int inNcat, int inbaseNcat, bool incondition_aware,
        bool inpolymorphism_aware, unsigned inprecision, bool inarithmetic,
        bool inmove_root_pop_size, bool inclamp_pop_sizes, bool inclamp_nuc_matrix,
        bool inclamp_corr_matrix, std::string infossilsfile, int inprior_cov_df, bool inuniq_kappa)
        : datafile(std::move(indatafile)),
          treefile(std::move(intreefile)),
          traitsfile(std::move(intraitsfile)),
          profilesfile(std::move(inprofiles)),
          fossilsfile(std::move(infossilsfile)),
          condition_aware(incondition_aware),
          polymorphism_aware(inpolymorphism_aware),
          precision(inprecision),
          arithmetic(inarithmetic),
          move_root_pop_size(inmove_root_pop_size),
          clamp_pop_sizes(inclamp_pop_sizes),
          clamp_nuc_matrix(inclamp_nuc_matrix),
          clamp_corr_matrix(inclamp_corr_matrix),
          prior_cov_df{inprior_cov_df},
          uniq_kappa{inuniq_kappa},
          baseNcat(inbaseNcat),
          Ncat(inNcat) {
        if (arithmetic) { std::cout << "Arithmetic mean instead of geodesic." << std::endl; }
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        if (PolymorphismAware()) { polydata = new PolyData(codondata, datafile); }

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        std::tuple<std::vector<std::vector<double>>, std::vector<size_t>> prefs{};
        if (Ncat <= 0) { Ncat = Nsite; }
        if (Ncat > Nsite) { Ncat = Nsite; }
        if (Ncat > 100) { Ncat = 100; }
        if (profilesfile != "Null") {
            prefs = open_preferences(profilesfile);
            clamp_profiles = true;
            move_root_pop_size = true;
            Ncat = static_cast<int>(std::get<0>(prefs).size());
            assert(static_cast<int>(std::get<1>(prefs).size()) == Nsite);
        }
        assert(Ncat <= Nsite);

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        std::cerr << "ncat : " << Ncat << '\n';
        std::cerr << "basencat : " << baseNcat << '\n';

        taxonset = codondata->GetTaxonSet();

        if (traitsfile != "Null") {
            taxon_traits = new TaxonTraits(traitsfile, *taxonset, polymorphism_aware);
        }

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        Nbranch = tree->nb_branches();
        branchdnds = new SimpleBranchArray<double>(*tree, 0.0);
        // Node ages
        nodeages = new NodeAges(*tree, fossilsfile);
        // Chronogram (diff between node ages)
        chronogram = new Chronogram(*nodeages);

        dimensions = 2;
        if (PolymorphismAware()) {
            dimensions++;
            assert(taxon_traits != nullptr);
        }
        if (taxon_traits != nullptr) { dimensions += taxon_traits->GetDim(); }
        prior_matrix = new PriorCovariance(dimensions, prior_cov_df, uniq_kappa);
        precision_matrix = new PrecisionMatrix(*prior_matrix);

        node_multivariate = new NodeMultivariateProcess(*chronogram, *precision_matrix, dimensions);

        // Branch Population size (brownian process)
        nodepopsize = new NodeProcess(*node_multivariate, dim_pop_size);
        branchpopsize = new BranchProcess(*nodepopsize, arithmetic);

        // Branch mutation rates (nbr of mutations per generation)
        nodemutrates = new NodeProcess(*node_multivariate, dim_mut_rate);
        nodemutrates->SlidingMove(-1);
        branchmutrates = new BranchProcess(*nodemutrates, arithmetic);

        if (PolymorphismAware()) {
            // Branch generation rate (nbr of generations per time)
            nodegentimes = new NodeProcess(*node_multivariate, dim_gen_time);
            branchgentimes = new BranchProcess(*nodegentimes, arithmetic);

            // Branch lengths (product of chronogram and mutation rate, divided by generation time)
            branchlength =
                new BranchwiseProductDivision(*chronogram, *branchmutrates, *branchgentimes);

            theta =
                new NodeProcessScaledMutationRate(theta_scale, nodepopsize, nodemutrates, Ntaxa);
        } else {
            // Branch lengths (product of branch mutation rate and chronogram)
            branchlength = new BranchwiseProduct(*chronogram, *branchmutrates);
        }

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix
        nucrelrate.assign(Nrr, 1.0 / Nrr);
        nucstat.assign(Nnuc, 1.0 / Nnuc);
        if (!clamp_nuc_matrix) {
            Random::DirichletSample(nucrelrate, std::vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
            Random::DirichletSample(nucstat, std::vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));
        }
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

        // selector, specifying which aa fitness array should be used for each site
        siteaafitnessarray =
            new MixtureSelector<std::vector<double>>(componentaafitnessarray, sitealloc);

        if (clamp_profiles) {
            for (int cat = 0; cat < Ncat; cat++) {
                (*componentaafitnessarray)[cat] = std::get<0>(prefs)[cat];
            }
            for (int site = 0; site < Nsite; site++) {
                (*sitealloc)[site] = static_cast<unsigned>(std::get<1>(prefs)[site]);
                assert(siteaafitnessarray->GetVal(site).size() == 20);
            }
        }

        // occupancy suff stats of site allocations (for resampling weights)
        occupancy = new OccupancySuffStat(Ncat);

        // codon matrices per branch and per site
        branchcomponentcodonmatrixarray = new MutSelNeCodonMatrixBidimArray(
            GetCodonStateSpace(), nucmatrix, componentaafitnessarray, branchpopsize->GetArray());

        // sub matrices per branch and per site
        branchsitecodonmatrixarray = new BranchComponentMatrixSelector<SubMatrix>(
            branchcomponentcodonmatrixarray, sitealloc, *tree);

        rootcomponentcodonmatrixarray = new AAMutSelNeCodonSubMatrixArray(
            GetCodonStateSpace(), nucmatrix, componentaafitnessarray, root_popsize);

        // sub matrices for root, across sites
        rootsitecodonmatrixarray =
            new RootComponentMatrixSelector<SubMatrix>(rootcomponentcodonmatrixarray, sitealloc);

        // The scaling factor of theta (4*Ne*u) since
        theta_scale = 1e-5;
        if (PolymorphismAware()) {
            poissonrandomfield = new PoissonRandomField(
                polydata->GetSampleSizeSet(), *GetCodonStateSpace(), precision);
            polyprocess = new PolyProcess(*GetCodonStateSpace(), *polydata, *poissonrandomfield,
                *siteaafitnessarray, *nucmatrix, *theta);
            taxoncomponentpolysuffstatbidimarray = new PolySuffStatBidimArray(Ntaxa, Ncat);
            taxonsitepolysuffstatbidimarray = new PolySuffStatBidimArray(Ntaxa, Nsite);
        }

        phyloprocess = new PhyloProcess(tree.get(), codondata, branchlength, nullptr,
            branchsitecodonmatrixarray, rootsitecodonmatrixarray, polyprocess);
        phyloprocess->Unfold();

        if (PolymorphismAware()) { theta->SetTaxonMap(&phyloprocess->GetTaxonMap()); }
        if (taxon_traits != nullptr) {
            node_multivariate->ClampLeaves(*taxon_traits, phyloprocess->GetTaxonMap());
        }
        branchcomponentpathsuffstatbidimarray = new PathSuffStatBidimArray(Nbranch, Ncat);
        branchsitepathsuffstatbidimarray = new PathSuffStatBidimArray(Nbranch, Nsite);
        rootcomponentpathsuffstatarray = new PathSuffStatArray(Ncat);
        rootsitepathsuffstatarray = new PathSuffStatArray(Nsite);

        scattersuffstat = new ScatterSuffStat(*tree);
        branchlengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
    }

    ~DatedNodeMutSelModel() override = default;

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        model_node(info, "nodeages", *nodeages);
        model_node(info, "root_popsize", root_popsize);
        model_node(info, "node_multivariate", *node_multivariate);
        model_node(info, "prior_cov_matrix", *prior_matrix);
        model_node(info, "precision_matrix", *precision_matrix);
        model_node(info, "nucrelrate", nucrelrate);
        model_node(info, "nucstat", nucstat);
        model_node(info, "basekappa", basekappa);
        model_node(info, "baseweight", *baseweight);
        model_node(info, "componentalloc", *componentalloc);
        model_node(info, "*basecenterarray", *basecenterarray);
        model_node(info, "*baseconcentrationarray", *baseconcentrationarray);
        model_node(info, "kappa", kappa);
        model_node(info, "weight", *weight);
        model_node(info, "componentaafitnessarray", *componentaafitnessarray);
        model_node(info, "sitealloc", *sitealloc);
        if (PolymorphismAware()) { model_node(info, "theta_scale", theta_scale); }
        model_node(info, "chronomove", chronomove);
        model_node(info, "chronoresamblesub", chronoresamblesub);

        // Descriptive statistics - Likelihood
        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        // Descriptive statistics - tree length
        model_stat(info, "length", [this]() { return 3 * branchlength->GetSum(); });
        // Descriptive statistics - dNdS predicted at mutation-selection balance
        model_stat(info, "PredictedDNDS", [this]() { return GetPredictedDNDS(); });
        // Descriptive statistics - Amino-acids preferences
        model_stat(info, "ncluster", [this]() { return GetNcluster(); });
        model_stat(info, "kappa", kappa);
        if (baseNcat > 1) {
            model_stat(info, "basencluster", [this]() { return GetBaseNcluster(); });
            model_stat(info, "basekappa", basekappa);
        }
        model_stat(info, "aaent", [this]() { return GetMeanAAEntropy(); });
        model_stat(info, "meanaaconc", [this]() { return GetMeanComponentAAConcentration(); });
        model_stat(info, "aacenterent", [this]() { return GetMeanComponentAAEntropy(); });
        // Descriptive statistics - Nucleotide matrix
        model_stat(info, "statent", [&]() { return Random::GetEntropy(nucstat); });
        model_stat(info, "rrent", [&]() { return Random::GetEntropy(nucrelrate); });
        // Descriptive statistics - Generator of the multivariate processes
        for (int i = 0; i < dimensions; i++) {
            model_stat(info, "PriorCovariance_" + std::to_string(i), prior_matrix->coeffRef(i));
            for (int j = 0; j <= i; j++) {
                model_stat(info, "Precision_" + std::to_string(i) + "_" + std::to_string(j),
                    precision_matrix->coeffRef(i, j));
            }
        }
        // Descriptive statistics - mean and var of mutation rate per generation
        model_stat(info, "BranchMutRatesMean", [this]() { return branchmutrates->GetMean(); });
        model_stat(info, "BranchMutRatesVar", [this]() { return branchmutrates->GetVar(); });
        // Descriptive statistics - mean and var of branch length
        model_stat(info, "BranchLengthMean", [this]() { return branchlength->GetMean(); });
        model_stat(info, "BranchLengthVar", [this]() { return branchlength->GetVar(); });
        // Descriptive statistics - mean and var of Ne
        model_stat(info, "BranchPopSizeMean", [this]() { return branchpopsize->GetMean(); });
        model_stat(info, "BranchPopSizeVar", [this]() { return branchpopsize->GetVar(); });
        if (PolymorphismAware()) {
            // Descriptive statistics - mean and var of the time per generation
            model_stat(info, "BranchGenTimesMean", [this]() { return branchgentimes->GetMean(); });
            model_stat(info, "BranchGenTimesVar", [this]() { return branchgentimes->GetVar(); });

            // Descriptive statistics - 4*Ne*u at the leaves of the tree
            model_stat(info, "ThetaScale", theta_scale);
            for (int taxon = 0; taxon < Ntaxa; taxon++) {
                model_stat(info, "*Theta_" + taxonset->GetTaxon(taxon), (*theta)[taxon]);
            }
        }
        // Descriptive statistics - mutation rate per generation at the root of the tree
        model_stat(info, "RootMutRate", (*nodemutrates)[tree->root()]);
        // Descriptive statistics - Ne under which the root sequence equilibrium is computed
        model_stat(info, "AncestralRootPopSize", root_popsize);
        // Descriptive statistics - Ne at the root of the tree (starting the multivariate process)
        model_stat(info, "RootPopSize", (*nodepopsize)[tree->root()]);
        // Descriptive statistics - time per generation at the root of the tree
        if (PolymorphismAware()) { model_stat(info, "RootGenTime", (*nodegentimes)[tree->root()]); }

        // Descriptive statistics - for each branch of the tree
        for (Tree::BranchIndex branch = 0; branch < tree->nb_branches(); branch++) {
            string b_name = tree->node_name(tree->node_index(branch));
            model_stat(info, "*BranchTime_" + b_name, (*chronogram)[branch]);
            model_stat(info, "*BranchMutRate_" + b_name, (*branchmutrates)[branch]);
            if (PolymorphismAware()) {
                model_stat(info, "*BranchGenTime_" + b_name, (*branchgentimes)[branch]);
            }
            model_stat(info, "*BranchLength_" + b_name, (*branchlength)[branch]);
            model_stat(info, "*BranchPopSize_" + b_name, (*branchpopsize)[branch]);
            model_stat(info, "*BranchdNdS_" + b_name, (*branchdnds)[branch]);
        }

        // Descriptive statistics - for chrono
        model_stat(
            info, "ChronoMove", [&]() { return chronomove / (chronomove + chronoresamblesub); });
        model_stat(info, "ChronoResampleSub",
            [&]() { return chronoresamblesub / (chronomove + chronoresamblesub); });
        for (auto const &chrono : moves_chrono_time) {
            model_stat(info, "Chrono" + chrono.first, moves_chrono_time[chrono.first]);
        }

        // Descriptive statistics - for each move acceptance
        for (auto const &mv : moves_tot) {
            model_stat(info, "Moves" + mv.first, moves_accepted[mv.first]);
        }
    }

    void ChronoStart(std::string const &name) {
        if (moves_chrono.count(name) == 0) {
            assert(moves_chrono_time.count(name) == 0);
            moves_chrono[name] = Chrono();
            moves_chrono_time[name] = 0.0;
        }
        moves_chrono[name].Reset();
        moves_chrono[name].Start();
    }

    void ChronoStop(std::string const &name) {
        moves_chrono[name].Stop();
        moves_chrono_time[name] += moves_chrono[name].GetTime();
    }

    std::string MoveName(std::string const &base, double tuning) {
        std::string name =
            base + "_" + std::to_string(tuning).substr(0, std::to_string(tuning).find('.') + 3);
        if (moves_accepted.count(name) == 0) {
            assert(moves_tot.count(name) == 0);
            moves_accepted[name] = 0.0;
            moves_tot[name] = 0.0;
        }
        return name;
    }

    void MoveAcceptation(std::string const &name, double accept_rate) {
        moves_tot[name]++;
        moves_accepted[name] += accept_rate;
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

    //! return tree
    const Tree &GetTree() const { return *tree; }

    //! return time of branch
    double GetBranchTime(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return chronogram->GetVal(tree->branch_index(node));
    };

    //! return branch length
    double GetBranchLength(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return branchlength->GetVal(tree->branch_index(node));
    };

    //! return contrast for a given node and dimension
    double GetContrast(Tree::NodeIndex node, int dim) const {
        assert(!tree->is_root(node));
        return node_multivariate->GetContrast(node)(dim);
    };

    //! return true if the model handles polymorphism
    bool PolymorphismAware() const { return polymorphism_aware; };

    //! return theta=4*Ne*u for a given node
    double GetTheta(Tree::NodeIndex node) const { return theta->GetNodeTheta(node); };

    //! return the value of the multivariate brownian process for a given node and a given dimension
    //! of the process
    double GetBrownianEntry(Tree::NodeIndex node, int dim) const {
        return node_multivariate->GetVal(node)(dim);
    }

    //! return number of dimensions of the multivariate brownian process
    int GetDimension() const { return dimensions; }

    std::string GetDimensionName(int dim) const {
        if (dim == dim_pop_size) {
            return "PopulationSize";
        } else if (dim == dim_mut_rate) {
            if (PolymorphismAware()) {
                return "MutationRatePerGeneration";
            } else {
                return "MutationRatePerTime";
            }
        } else if (dim == dim_gen_time and PolymorphismAware()) {
            return "GenerationTime";
        } else {
            assert(taxon_traits != nullptr);
            return "Traits" + taxon_traits->GetHeader(dim);
        }
    }
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
        for (Tree::BranchIndex branch = 0; branch < tree->nb_branches(); branch++) {
            branchcomponentcodonmatrixarray->UpdateRowNe(branch, branchpopsize->GetVal(branch));
        }
        rootcomponentcodonmatrixarray->UpdateNe(root_popsize);
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatricesNoFitnessRecomput() {
        branchcomponentcodonmatrixarray->UpdateCodonMatricesNoFitnessRecomput();
        rootcomponentcodonmatrixarray->UpdateCodonMatricesNoFitnessRecomput();
    }

    //! \brief tell codon matrices for cat i and across branch that their parameters have
    //! changed and that they should be updated
    void UpdateCatCodonMatrix(int i) {
        branchcomponentcodonmatrixarray->UpdateColCodonMatrices(i);
        rootcomponentcodonmatrixarray->UpdateCodonMatrices(i);
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatricesNoFitnessRecomput() {
        UpdateNucMatrix();
        UpdateCodonMatricesNoFitnessRecomput();
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
        branchmutrates->Update();
        if (PolymorphismAware()) { branchgentimes->Update(); }
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
        branchmutrates->UpdateLocal(node);
        if (PolymorphismAware()) { branchgentimes->UpdateLocal(node); }
        branchlength->UpdateLocal(node);
    }

    //! \brief Update the branch Ne around the focal node.
    //!
    //! Update needed when the Ne (NodeProcess) of the focal node is changed.
    void UpdateLocalBranchPopSize(Tree::NodeIndex node) {
        assert(!tree->is_root(node) or clamp_profiles);
        branchpopsize->UpdateLocal(node);
        UpdateBranchPopSize(node);
        for (Tree::NodeIndex const &child : tree->children(node)) { UpdateBranchPopSize(child); }
    }

    void UpdateBranchPopSize(Tree::NodeIndex node) {
        Tree::BranchIndex branch = tree->branch_index(node);
        branchcomponentcodonmatrixarray->UpdateRowNe(branch, branchpopsize->GetVal(branch));
    }

    //! \brief Update the root Ne.
    //!
    //! Update needed when the Ne of the root is changed.
    void UpdateRootPopSize() { rootcomponentcodonmatrixarray->UpdateNe(root_popsize); }

    void Update() {
        UpdateBranches();
        UpdateBaseOccupancies();
        UpdateOccupancies();
        UpdateMatrices();
        UpdateStats();
    }

    void UpdateStats() {
        if (PolymorphismAware()) { theta->Update(); }
        for (Tree::BranchIndex b = 0; b < Nbranch; b++) { (*branchdnds)[b] = GetPredictedDNDS(b); }
    }

    void PostPred(std::string name) {
        Update();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior (up to some constant)
    double GetLogPrior() const {
        double total = PrecisionMatrixLogProb();
        total += NodeMultivariateLogPrior();
        total += RootMultivariateLogPrior();

        if (!clamp_nuc_matrix) { total += NucRatesLogPrior(); }

        if (!clamp_profiles) {
            if (baseNcat > 1) {
                total += BaseStickBreakingHyperLogPrior();
                total += BaseStickBreakingLogPrior();
            }
            total += BaseLogPrior();
            total += StickBreakingHyperLogPrior();
            total += StickBreakingLogPrior();
            total += AALogPrior();
        }

        if (PolymorphismAware()) { total += ThetaLogPrior(); }
        return total;
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    //! log prob of precision matrix
    double PrecisionMatrixLogProb() const {
        return prior_matrix->GetLogProb() + precision_matrix->GetLogProb(*prior_matrix);
    }

    // Multivariate prior
    //! log prior over branch rates (brownian process)
    double NodeMultivariateLogPrior() const { return node_multivariate->GetLogProb(); }

    //! log prior of branch rate (brownian process) around of focal node
    double LocalNodeMultivariateLogPrior(Tree::NodeIndex node) const {
        return node_multivariate->GetLocalLogProb(node);
    }

    //! log prior of
    double RootMultivariateLogPrior() const { return node_multivariate->GetLogProb(tree->root()); }

    //! log prior over theta
    double ThetaLogPrior() const { return 0.0; }

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
        assert(!std::isinf(total));
        return total;
    }

    //! log prior over base center and concentration parameters of component k of
    //! base distribution
    double BaseLogPrior(int k) const {
        double total = 0;
        total += basecenterarray->GetLogProb(k);
        total += baseconcentrationarray->GetLogProb(k);
        assert(!std::isinf(total));
        return total;
    }

    //! log prior of amino-acid fitness profiles
    double AALogPrior() const { return componentaafitnessarray->GetLogProb(); }

    //! log prior of amino-acid fitness profile k
    double AALogPrior(int k) const { return componentaafitnessarray->GetLogProb(k); }


    //-------------------
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics for substitution mappings across sites
    void CollectSitePathSuffStat() {
        branchsitepathsuffstatbidimarray->Clear();
        rootsitepathsuffstatarray->Clear();
        branchsitepathsuffstatbidimarray->AddSuffStat(*phyloprocess, *rootsitepathsuffstatarray);
    }

    //! gather site-specific sufficient statistics component-wise
    void CollectComponentPathSuffStat() {
        branchcomponentpathsuffstatbidimarray->Clear();
        branchcomponentpathsuffstatbidimarray->Add(*branchsitepathsuffstatbidimarray, *sitealloc);
        rootcomponentpathsuffstatarray->Clear();
        rootcomponentpathsuffstatarray->Add(*rootsitepathsuffstatarray, *sitealloc);
    }

    //! collect sufficient statistics at the tips of the tree
    void CollectSitePolySuffStat() {
        if (PolymorphismAware()) {
            taxonsitepolysuffstatbidimarray->Clear();
            phyloprocess->AddPolySuffStat(*taxonsitepolysuffstatbidimarray);
        }
    }

    //! gather site-specific tips of the tree sufficient statistics component-wise
    void CollectComponentPolySuffStat() {
        if (PolymorphismAware()) {
            taxoncomponentpolysuffstatbidimarray->Clear();
            taxoncomponentpolysuffstatbidimarray->Add(*taxonsitepolysuffstatbidimarray, *sitealloc);
        }
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        branchlengthpathsuffstatarray->Clear();
        branchlengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! collect sufficient statistics for sampling the correlation matrix of the brownian process
    void CollectScatterSuffStat() {
        scattersuffstat->Clear();
        scattersuffstat->AddSuffStat(*node_multivariate);
    }

    //! collect suff stats for moving center and concentration parameters of the
    //! base mixture
    void CollectBaseSuffStat() {
        basesuffstatarray->Clear();
        componentaafitnessarray->AddSuffStat(*basesuffstatarray, *componentalloc);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! return log prob only at the tips due to polymorphism of the substitution mapping,
    //! as a function of the current codon substitution process
    double PolySuffStatLogProb() const {
        //! sum over all components to get log prob
        if (PolymorphismAware()) {
            return taxoncomponentpolysuffstatbidimarray->GetLogProb(
                *poissonrandomfield, *componentaafitnessarray, *nucmatrix, *theta);
        } else {
            return 0;
        }
    }

    //! return log prob only at the tips due to polymorphism of the substitution
    //! mapping, over sites allocated to component k of the mixture
    double ComponentPolySuffStatLogProb(int k) const {
        // sum over all sites allocated to component k
        if (PolymorphismAware()) {
            double tot = 0;
            for (int taxon = 0; taxon < Ntaxa; taxon++) {
                tot += taxoncomponentpolysuffstatbidimarray->GetVal(taxon, k).GetLogProb(
                    *poissonrandomfield, componentaafitnessarray->GetVal(k), *nucmatrix,
                    theta->GetTheta(taxon));
            }
            return tot;
        } else {
            return 0.0;
        }
    }

    //! return log prob only at the tips due to polymorphism of the substitution
    //! mapping, over sites allocated to component k of the mixture
    double TaxonPolySuffStatLogProb(int taxon) const {
        // sum over all sites allocated to component k
        if (PolymorphismAware()) {
            double tot = 0;
            double d_theta = theta->GetTheta(taxon);
            for (int k = 0; k < Ncat; k++) {
                tot += taxoncomponentpolysuffstatbidimarray->GetVal(taxon, k).GetLogProb(
                    *poissonrandomfield, componentaafitnessarray->GetVal(k), *nucmatrix, d_theta);
            }
            return tot;
        } else {
            return 0.0;
        }
    }

    //! return log prob only at the tips due to polymorphism of the substitution
    //! mapping, for a given sites if allocated to component cat of the mixture
    double SitePolySuffStatLogProbGivenComponent(int site, int cat) const {
        if (PolymorphismAware()) {
            double tot = 0;
            for (int taxon = 0; taxon < Ntaxa; taxon++) {
                tot += taxonsitepolysuffstatbidimarray->GetVal(taxon, site)
                           .GetLogProb(*poissonrandomfield, componentaafitnessarray->GetVal(cat),
                               *nucmatrix, theta->GetTheta(taxon));
            }
            return tot;
        } else {
            return 0.0;
        }
    }

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return branchcomponentpathsuffstatbidimarray->GetLogProb(*branchcomponentcodonmatrixarray) +
               RootPathSuffStatLogProb() + PolySuffStatLogProb();
    }

    double RootPathSuffStatLogProb() const {
        return rootcomponentpathsuffstatarray->GetLogProb(*rootcomponentcodonmatrixarray);
    }

    //! return log prob of the substitution mappings, for a given sites if allocated to component
    //! cat of the mixture
    double SitePathSuffStatLogProbGivenComponent(int site, int cat) const {
        return branchsitepathsuffstatbidimarray->GetColLogProb(
                   site, *branchcomponentcodonmatrixarray, cat) +
               rootsitepathsuffstatarray->GetVal(site).GetLogProb(
                   rootcomponentcodonmatrixarray->GetVal(cat)) +
               SitePolySuffStatLogProbGivenComponent(site, cat);
    }

    //! return log prob of the substitution mappings over sites allocated to
    //! component k of the mixture
    double ComponentPathSuffStatLogProb(int cat) const {
        return branchcomponentpathsuffstatbidimarray->GetColLogProb(
                   cat, *branchcomponentcodonmatrixarray) +
               rootcomponentpathsuffstatarray->GetVal(cat).GetLogProb(
                   rootcomponentcodonmatrixarray->GetVal(cat)) +
               ComponentPolySuffStatLogProb(cat);
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
        if (GetTree().is_root(node)) {
            return NodeMultivariateLogPrior() + BranchLengthSuffStatLogProb();
        } else {
            return LocalNodeMultivariateLogPrior(node) + LocalBranchLengthSuffStatLogProb(node);
        }
    }

    //! \brief log prob to be recomputed when moving branch rates (brownian process) around of focal
    //! node
    double LocalNodeRateLogProb(Tree::NodeIndex node, bool mut_rate_change) const {
        double tot = LocalNodeMultivariateLogPrior(node) + LocalBranchLengthSuffStatLogProb(node);
        if (mut_rate_change and PolymorphismAware() and tree->is_leaf(node)) {
            tot += TaxonPolySuffStatLogProb(phyloprocess->GetTaxonMap().NodeToTaxon(node));
        }
        return tot;
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
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! the length of a given branch
    double BranchLengthSuffStatLogProb(Tree::BranchIndex branch) const {
        return branchlengthpathsuffstatarray->GetVal(branch).GetLogProb(
            branchlength->GetVal(branch));
    }

    //! \brief return log prob of current substitution mapping (on all branches), as a function of
    //! the length all branches
    double BranchLengthSuffStatLogProb() const {
        return branchlengthpathsuffstatarray->GetLogProb(*branchlength);
    }

    //! \brief log prob to be recomputed when moving Ne (brownian process) around of focal node
    double LocalNodePopSizeLogProb(Tree::NodeIndex node) const {
        return LocalNodeMultivariateLogPrior(node) + LocalNodePopSizeSuffStatLogProb(node);
    }

    //! \brief log prob factor (without prior) to be recomputed when moving Ne (brownian process)
    //! around of focal node
    double LocalNodePopSizeSuffStatLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        // for all children
        for (auto const &child : tree->children(node)) {
            tot += NodePopSizeSuffStatLogProb(tree->branch_index(child));
        }
        if (tree->is_root(node)) {
            // for the root we use the rootpathsuffstat
            tot += RootPathSuffStatLogProb();
        } else {
            // for the branch attached to the node
            tot += NodePopSizeSuffStatLogProb(tree->branch_index(node));
        }
        if (PolymorphismAware() and tree->is_leaf(node)) {
            tot += TaxonPolySuffStatLogProb(phyloprocess->GetTaxonMap().NodeToTaxon(node));
        }
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! Ne of a given branch
    double NodePopSizeSuffStatLogProb(Tree::BranchIndex branch) const {
        return branchcomponentpathsuffstatbidimarray->GetRowLogProb(
            branch, *branchcomponentcodonmatrixarray);
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
    void Move() {
        for (auto &mv : moves_tot) {
            moves_accepted[mv.first] = 0.0;
            moves_tot[mv.first] = 0.0;
        }
        for (auto &chrono_time : moves_chrono_time) { chrono_time.second = 0.0; }

        Chrono chrono;
        chrono.Start();
        ResampleSub(1.0);
        chrono.Stop();
        chronoresamblesub = chrono.GetTime();

        chrono.Reset();
        chrono.Start();
        MoveParameters(20);
        chrono.Stop();
        chronomove = chrono.GetTime();

        for (auto &mv : moves_accepted) { moves_accepted[mv.first] /= moves_tot[mv.first]; }
        for (auto &chrono_time : moves_chrono_time) { chrono_time.second /= chronomove; }
        UpdateStats();
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        UpdateMatrices();
        phyloprocess->Move(frac);
        assert(CheckMapping());
    }

    bool CheckMapping() const {
        for (int taxon = 0; taxon < Ntaxa; taxon++) {
            for (int site = 0; site < Nsite; site++) {
                int sub_state = phyloprocess->GetPathState(taxon, site);
                int data_state = codondata->GetState(taxon, site);
                if (data_state == -1) { continue; }
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
                    return false;
                }
            }
        }
        return true;
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        CollectSitePolySuffStat();
        for (int rep = 0; rep < nrep; rep++) {
            CollectComponentPolySuffStat();
            CollectLengthSuffStat();

            ChronoStart("NodeAges");
            MoveNodeAges(0.1, 3, true);
            MoveNodeAges(0.1, 3, false);
            MoveNodeAges(0.02, 3, true);
            MoveNodeAges(0.02, 3, false);
            ChronoStop("NodeAges");

            ChronoStart("NodeRates");
            MoveNodeRates(0.4, 3, true);
            MoveNodeRates(0.4, 3, false);
            MoveNodeRates(0.2, 3, true);
            MoveNodeRates(0.2, 3, false);
            ChronoStop("NodeRates");

            CollectSitePathSuffStat();
            CollectComponentPathSuffStat();
            if (!clamp_nuc_matrix) {
                ChronoStart("NucRates");
                MoveNucRates();
                ChronoStop("NucRates");
            }
            if (!clamp_pop_sizes) {
                if (move_root_pop_size) {
                    ChronoStart("NodeRootPopSizes");
                    MoveRootPopSize(0.4, 6);
                    MoveRootPopSize(0.2, 6);
                    ChronoStop("NodeRootPopSizes");
                }
                ChronoStart("NodePopSizes");
                MoveNodePopSizes(0.4, 3, true);
                MoveNodePopSizes(0.4, 3, false);
                MoveNodePopSizes(0.2, 3, true);
                MoveNodePopSizes(0.2, 3, false);
                ChronoStop("NodePopSizes");
            }

            ChronoStart("NodeTraits");
            MoveNodeTraits(0.5, 3, true);
            MoveNodeTraits(0.5, 3, false);
            MoveNodeTraits(0.05, 3, true);
            MoveNodeTraits(0.05, 3, false);
            ChronoStop("NodeTraits");

            if (!clamp_corr_matrix) {
                ChronoStart("PrecisionMatrix");
                CollectScatterSuffStat();
                SamplePrecisionMatrix();
                MovePriorMatrix(0.1, 3);
                MovePriorMatrix(0.01, 3);
                ChronoStop("PrecisionMatrix");
            }

            if (PolymorphismAware()) {
                ChronoStart("Theta");
                MoveTheta();
                ChronoStop("Theta");
            }

            if (!clamp_profiles) {
                ChronoStart("AAMixture");
                MoveAAMixture(2);
                ChronoStop("AAMixture");
                ChronoStart("Base");
                MoveBase(2);
                ChronoStop("Base");
            }
        }
    }

    void SamplePrecisionMatrix() {
        scattersuffstat->SamplePrecisionMatrix(*precision_matrix, *prior_matrix);
    };

    //! MH moves on the invert wishart matrix (prior of the covariance matrix)
    void MovePriorMatrix(double tuning, int nrep) {
        for (int i = 0; i < (uniq_kappa ? 1 : dimensions); ++i) {
            Move::Scaling((*prior_matrix)(i), tuning, nrep,
                &DatedNodeMutSelModel::PrecisionMatrixLogProb, &DatedNodeMutSelModel::NoUpdate,
                this);
        }
    }

    //! MH moves on branch ages
    void MoveNodeAges(double tuning, int nrep, bool leaves_to_root) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node :
                leaves_to_root ? tree->leaves_root_to_iter() : tree->root_to_leaves_iter()) {
                if (!tree->is_leaf(node)) { MoveNodeAge(node, tuning); }
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
        MoveAcceptation(MoveName("NodeAge", tuning), accept);
    }

    //! MH moves on traits (brownian process)
    void MoveNodeTraits(double tuning, int nrep, bool leaves_to_root) {
        if (taxon_traits == nullptr or taxon_traits->GetDim() == 0) { return; }
        for (int rep = 0; rep < nrep; rep++) {
            for (int trait_dim = 0; trait_dim < taxon_traits->GetDim(); trait_dim++) {
                for (Tree::NodeIndex node :
                    leaves_to_root ? tree->leaves_root_to_iter() : tree->root_to_leaves_iter()) {
                    if (!tree->is_leaf(node) or
                        !taxon_traits->DataPresence(
                            phyloprocess->GetTaxonMap().NodeToTaxon(node), trait_dim)) {
                        MoveNodeTrait(node, tuning, trait_dim);
                    }
                }
            }
        }
    }

    //! MH moves on traits (brownian process) for a focal node
    void MoveNodeTrait(Tree::NodeIndex node, double tuning, int trait_dim) {
        int multi_dim = taxon_traits->TraitDimToMultivariateDim(trait_dim);
        double logratio = -LocalNodeMultivariateLogPrior(node);

        double m = tuning * (Random::Uniform() - 0.5);
        (*node_multivariate)[node](multi_dim) += m;

        logratio += LocalNodeMultivariateLogPrior(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) { (*node_multivariate)[node](multi_dim) -= m; }
        MoveAcceptation(MoveName("NodeTrait", tuning), accept);
    }

    //! MH moves on branch rates (brownian process)
    void MoveNodeRates(double tuning, int nrep, bool leaves_to_root) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node :
                leaves_to_root ? tree->leaves_root_to_iter() : tree->root_to_leaves_iter()) {
                MoveNodeRate(node, tuning, true);
                if (PolymorphismAware() and
                    (!tree->is_leaf(node) or !taxon_traits->GenTimePresence(
                                                 phyloprocess->GetTaxonMap().NodeToTaxon(node)))) {
                    MoveNodeRate(node, tuning, false);
                }
            }
        }
    }

    //! MH moves on branch rates (brownian process) for a focal node
    void MoveNodeRate(Tree::NodeIndex node, double tuning, bool mut_rates) {
        double logratio = -LocalNodeRateLogProb(node, mut_rates);

        double m = tuning * (Random::Uniform() - 0.5);
        (mut_rates ? nodemutrates : nodegentimes)->SlidingMove(node, m);
        UpdateLocalBranchRates(node);

        logratio += LocalNodeRateLogProb(node, mut_rates);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            (mut_rates ? nodemutrates : nodegentimes)->SlidingMove(node, -m);
            UpdateLocalBranchRates(node);
        }
        MoveAcceptation(MoveName("NodeRate", tuning), accept);
    }


    //! MH moves on branch Ne (brownian process)
    void MoveNodePopSizes(double tuning, int nrep, bool leaves_to_root) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node :
                leaves_to_root ? tree->leaves_root_to_iter() : tree->root_to_leaves_iter()) {
                if (!tree->is_root(node) or clamp_profiles) { MoveNodePopSize(node, tuning); }
            }
        }
    }

    //! MH moves on branch Ne (brownian process) for a focal node
    void MoveNodePopSize(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodePopSizeLogProb(node);

        double m = tuning * (Random::Uniform() - 0.5);
        nodepopsize->SlidingMove(node, m);
        UpdateLocalBranchPopSize(node);

        logratio += LocalNodePopSizeLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            nodepopsize->SlidingMove(node, -m);
            UpdateLocalBranchPopSize(node);
        }
        MoveAcceptation(MoveName("NodePopSize", tuning), accept);
    }

    //! MH moves on branch Ne (brownian process)
    void MoveRootPopSize(double tuning, int nrep) {
        double rate = Move::Scaling(root_popsize, tuning, nrep,
            &DatedNodeMutSelModel::RootPathSuffStatLogProb,
            &DatedNodeMutSelModel::UpdateRootPopSize, this);
        MoveAcceptation(MoveName("RootPopSize", tuning), rate);
    }

    //! MH move on base mixture
    void MoveBase(int nrep) {
        if (baseNcat > 1) { ResampleBaseAlloc(); }
        MoveBaseMixture(nrep);
    }


    //! MH move on theta
    void MoveTheta() {
        MoveTheta(1.0, 10);
        MoveTheta(0.3, 10);
    }

    //! MH move on theta
    void MoveTheta(double tuning, int nrep) {
        double rate = Move::Scaling(theta_scale, tuning, nrep, &DatedNodeMutSelModel::ThetaLogProb,
            &DatedNodeMutSelModel::NoUpdate, this);
        MoveAcceptation(MoveName("ThetaScale", tuning), rate);
    }

    //! MH move on nucleotide rate parameters
    void MoveNucRates(double tuning, int n, int nrep) {
        double rate =
            Move::Profile(nucrelrate, tuning, n, nrep, &DatedNodeMutSelModel::NucRatesLogProb,
                &DatedNodeMutSelModel::UpdateMatricesNoFitnessRecomput, this);
        MoveAcceptation(MoveName("NucRates_" + std::to_string(n), tuning), rate);
    }


    //! MH move on nucleotide rate parameters
    void MoveNucStat(double tuning, int n, int nrep) {
        double rate =
            Move::Profile(nucstat, tuning, n, nrep, &DatedNodeMutSelModel::NucRatesLogProb,
                &DatedNodeMutSelModel::UpdateMatricesNoFitnessRecomput, this);
        MoveAcceptation(MoveName("NucStat_" + std::to_string(n), tuning), rate);
    }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        MoveNucRates(0.1, 1, 3);
        MoveNucRates(0.03, 3, 3);
        MoveNucRates(0.01, 3, 3);
        MoveNucStat(0.1, 1, 3);
        MoveNucStat(0.01, 1, 3);
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
        branchcomponentcodonmatrixarray->UpdateColCodonMatrices(*occupancy);
        rootcomponentcodonmatrixarray->UpdateCodonMatrices(*occupancy);
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
        MoveAAGamma(0.4, nrep);
        MoveAAGamma(0.1, nrep);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    void MoveAA(double tuning, int n, int nrep) {
        auto name = MoveName("AA_" + std::to_string(n), tuning);

        double bk[Naa];
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                std::vector<double> &aa = (*componentaafitnessarray)[i];
                for (int rep = 0; rep < nrep; rep++) {
                    for (int l = 0; l < Naa; l++) { bk[l] = aa[l]; }
                    double deltalogprob = -AALogPrior(i) - ComponentPathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                    deltalogprob += loghastings;
                    UpdateCatCodonMatrix(i);
                    deltalogprob += AALogPrior(i) + ComponentPathSuffStatLogProb(i);
                    bool accepted = (log(Random::Uniform()) < deltalogprob);
                    if (!accepted) {
                        for (int l = 0; l < Naa; l++) { aa[l] = bk[l]; }
                        UpdateCatCodonMatrix(i);
                    }
                    MoveAcceptation(name, accepted);
                }
            }
        }
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
    void MoveAAGamma(double tuning, int nrep) {
        auto name = MoveName("AAGamma", tuning);

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

                    UpdateCatCodonMatrix(i);

                    deltalogprob +=
                        GammaAALogPrior(x, aacenter, aaconc) + ComponentPathSuffStatLogProb(i);

                    bool accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        bkaa = aa;
                        bkx = x;
                        bkz = z;
                    } else {
                        aa = bkaa;
                        x = bkx;
                        z = bkz;
                        UpdateCatCodonMatrix(i);
                    }
                    MoveAcceptation(name, accepted);
                }
            }
        }
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

        for (int cat = 0; cat < Ncat; cat++) {
            double tmp = SitePathSuffStatLogProbGivenComponent(site, cat);
            postprob[cat] = tmp;
            if ((!cat) || (max < tmp)) { max = tmp; }
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
    void MoveKappa(double tuning, int nrep) {
        double rate =
            Move::Scaling(kappa, tuning, nrep, &DatedNodeMutSelModel::StickBreakingHyperLogProb,
                &DatedNodeMutSelModel::NoUpdate, this);
        MoveAcceptation(MoveName("Kappa", tuning), rate);
    }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa() {
        MoveKappa(1.0, 10);
        MoveKappa(0.3, 10);
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
            MoveBaseCenters(0.3, 3);
            MoveBaseConcentrations(1.0);
            MoveBaseConcentrations(0.4);
        }
    }

    //! MCMC module for moving the center parameters of the components of the the
    //! base mixture
    void MoveBaseCenters(double tuning, int n) {
        auto name = MoveName("BaseCenters_" + std::to_string(n), tuning);

        std::vector<double> bk(Naa, 0);
        for (int k = 0; k < baseNcat; k++) {
            if (baseoccupancy->GetVal(k)) {
                std::vector<double> &aa = (*basecenterarray)[k];
                bk = aa;
                double deltalogprob = -BaseLogProb(k);
                double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                deltalogprob += loghastings;
                deltalogprob += BaseLogProb(k);
                bool accepted = (log(Random::Uniform()) < deltalogprob);
                if (!accepted) { aa = bk; }
                MoveAcceptation(name, accepted);
            }
        }
    }

    //! MCMC module for moving the concentration parameters of the components of
    //! the the base mixture
    void MoveBaseConcentrations(double tuning) {
        auto name = MoveName("BaseConc", tuning);

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
                bool accepted = (log(Random::Uniform()) < deltalogprob);
                if (!accepted) { c = bk; }
                MoveAcceptation(name, accepted);
            }
        }
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
    void MoveBaseKappa(double tuning, int nrep) {
        double rate = Move::Scaling(basekappa, tuning, nrep,
            &DatedNodeMutSelModel::BaseStickBreakingHyperLogProb, &DatedNodeMutSelModel::NoUpdate,
            this);
        MoveAcceptation(MoveName("BaseKappa", tuning), rate);
    }

    //! MH move on basekappa, concentration parameter of the base mixture
    void MoveBaseKappa() {
        MoveBaseKappa(2.0, 10);
        MoveBaseKappa(1.0, 10);
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
    double GetMeanAAEntropy() const {
        double aaent = 0;
        for (int k = 0; k < Ncat; k++) {
            if (occupancy->GetVal(k)) {
                aaent += occupancy->GetVal(k) * componentaafitnessarray->GetMeanEntropy(k);
            }
        }
        return aaent / GetNsite();
    }

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

    double GetPredictedDNDS(Tree::BranchIndex branch) const {
        double dn{0.}, dn0{0.};
        for (int k = 0; k < Ncat; k++) {
            if (occupancy->GetVal(k)) {
                double cat_dn{0}, cat_dn0{0};
                std::tie(cat_dn, cat_dn0) =
                    branchcomponentcodonmatrixarray->GetVal(branch, k).GetFlowDNDS();
                dn += occupancy->GetVal(k) * cat_dn;
                dn0 += occupancy->GetVal(k) * cat_dn0;
            }
        }
        return dn / dn0;
    }

    double GetPredictedDNDS() const { return branchdnds->GetMean(); }

    const std::vector<double> &GetProfile(int i) const { return siteaafitnessarray->GetVal(i); }

    void ToStream(std::ostream &os) { os << *this; }
};

std::istream &operator>>(std::istream &is, std::unique_ptr<DatedNodeMutSelModel> &m) {
    std::string model_name, datafile, treefile, traitsfile, profilesfile, fossilsfile;
    int Ncat, baseNcat, prior_cov_df;
    bool condition_aware, polymorphism_aware, arithmetic, move_root_pop_size, clamp_pop_sizes,
        clamp_nuc_matrix, clamp_corr_matrix, uniq_kappa;
    unsigned precision;

    is >> model_name;
    if (model_name != "DatedNodeMutSelModel") {
        std::cerr << "Expected DatedNodeMutSelModel for model name, got " << model_name << "\n";
        exit(1);
    }

    is >> datafile >> treefile >> traitsfile >> profilesfile >> fossilsfile;
    is >> Ncat >> baseNcat;
    is >> condition_aware >> polymorphism_aware >> precision >> arithmetic >> move_root_pop_size >>
        clamp_pop_sizes >> clamp_nuc_matrix >> clamp_corr_matrix >> prior_cov_df >> uniq_kappa;
    m = std::make_unique<DatedNodeMutSelModel>(datafile, treefile, traitsfile, profilesfile, Ncat,
        baseNcat, condition_aware, polymorphism_aware, precision, arithmetic, move_root_pop_size,
        clamp_pop_sizes, clamp_nuc_matrix, clamp_corr_matrix, fossilsfile, prior_cov_df,
        uniq_kappa);
    Tracer tracer{*m};
    tracer.read_line(is);
    m->Update();
    return is;
}

std::ostream &operator<<(std::ostream &os, DatedNodeMutSelModel &m) {
    Tracer tracer{m};
    os << "DatedNodeMutSelModel" << '\t';
    os << m.datafile << '\t';
    os << m.treefile << '\t';
    assert(!m.traitsfile.empty());
    os << m.traitsfile << '\t';
    assert(!m.profilesfile.empty());
    os << m.profilesfile << '\t';
    assert(!m.fossilsfile.empty());
    os << m.fossilsfile << '\t';
    os << m.Ncat << '\t';
    os << m.baseNcat << '\t';
    os << m.condition_aware << '\t';
    os << m.polymorphism_aware << '\t';
    os << m.precision << '\t';
    os << m.arithmetic << '\t';
    os << m.move_root_pop_size << '\t';
    os << m.clamp_pop_sizes << '\t';
    os << m.clamp_nuc_matrix << '\t';
    os << m.clamp_corr_matrix << '\t';
    os << m.prior_cov_df << '\t';
    os << m.uniq_kappa << '\t';
    tracer.write_line(os);
    return os;
}
