#pragma once

#include "Chronogram.hpp"
#include "Move.hpp"
#include "MultivariateProcess.hpp"
#include "ScatterSuffStat.hpp"
#include "TaxonTraits.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"

/**
 * \brief Traits are modeled as a brownian process along the tree (like CoEvol)
 *
 */


class DatedNodeModel : public ChainComponent {
    // tree and data
    std::string treefile, traitsfile;
    std::unique_ptr<const Tree> tree;
    TaxonSet *taxonset{nullptr};
    TaxonMap *taxonmap{nullptr};
    TaxonTraits *taxon_traits{nullptr};

    int Ntaxa;
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
    EMatrix *cov_matrix;
    NodeMultivariateProcess *node_multivariate;
    ScatterSuffStat *scattersuffstat;
    std::vector<double> root_trait;

  public:
    friend std::ostream &operator<<(std::ostream &os, DatedNodeModel &m);

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    DatedNodeModel(
        std::string intreefile, std::string &intraitsfile, int inprior_cov_df, bool inuniq_kappa)
        : treefile(std::move(intreefile)),
          traitsfile(std::move(intraitsfile)),
          prior_cov_df{inprior_cov_df},
          uniq_kappa{inuniq_kappa} {
        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        // Node ages
        nodeages = new NodeAges(*tree, parser.get_tree());
        chronogram = new Chronogram(*nodeages);
        chronogram->Scale();
        std::vector<std::string> leave_names;
        for (auto const &n : tree->root_to_leaves_iter()) {
            if (tree->is_leaf(n)) { leave_names.push_back(tree->node_name(n)); }
        }
        std::cout << "Tree has " << tree->nb_nodes() << " nodes and " << leave_names.size()
                  << " leaves\n";
        std::cout << "Total tree size: " << chronogram->GetSum() << std::endl;

        taxonset = new TaxonSet(leave_names);
        taxonmap = new TaxonMap(tree.get(), taxonset);
        Ntaxa = taxonmap->GetNtaxa();
        taxon_traits = new TaxonTraits(traitsfile, *taxonset, false, 0);

        // Precision matrix and prior
        dimensions = taxon_traits->GetDim();
        if (dimensions == 1) { uniq_kappa = true; }
        prior_matrix = new PriorCovariance(dimensions, prior_cov_df, uniq_kappa);
        precision_matrix = new PrecisionMatrix(*prior_matrix);
        cov_matrix = new EMatrix(dimensions, dimensions);
        precision_matrix->UpdateCovarianceMatrix(*cov_matrix);

        node_multivariate = new NodeMultivariateProcess(*chronogram, *precision_matrix, dimensions);
        if (taxon_traits != nullptr) { node_multivariate->ClampLeaves(*taxon_traits, *taxonmap); }
        // Suff Stats
        scattersuffstat = new ScatterSuffStat(*tree);
        root_trait = std::vector<double>(dimensions, 0);
    }

    virtual ~DatedNodeModel() = default;

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        model_node(info, "node_multivariate", *node_multivariate);
        model_node(info, "prior_cov_matrix", *prior_matrix);
        model_node(info, "precision_matrix", *precision_matrix);
        model_node(info, "cov_matrix", *cov_matrix);

        model_stat(info, "n_taxa", [&]() { return Ntaxa; });
        model_stat(info, "lnprob", [&]() { return DatedNodeModel::GetLogProb(); });
        model_stat(info, "priorLnprob", [&]() { return prior_matrix->GetLogProb(); });
        model_stat(
            info, "precisionLnprob", [&]() { return precision_matrix->GetLogProb(*prior_matrix); });
        model_stat(info, "multiLnprob", [&]() { return node_multivariate->GetLogProb(); });
        for (int i = 0; i < dimensions; i++) {
            model_stat(info, "PriorSigma_" + GetDimensionName(i), prior_matrix->coeffRef(i));
            // Trait at the root of the tree
            model_stat(info, "RootTrait_" + GetDimensionName(i), root_trait[i]);
            model_stat(info, "Var_" + GetDimensionName(i), cov_matrix->coeffRef(i, i));
            for (int j = 0; j <= i; j++) {
                model_stat(info, "Precision_" + GetDimensionName(i) + "_" + GetDimensionName(j),
                    precision_matrix->coeffRef(i, j));
                model_stat(info, "Cov_" + GetDimensionName(i) + "_" + GetDimensionName(j),
                    cov_matrix->coeffRef(i, j));
            }
        }
        model_stat(info, "ChronogramLength", [&]() { return GetChronogramLength(); });
    }

    //! return tree
    const Tree &GetTree() const { return *tree; }

    //! return time of branch
    double GetBranchTime(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return chronogram->GetVal(tree->branch_index(node));
    };

    //! return precision matrix
    PrecisionMatrix GetPrecisionMatrix() const { return *precision_matrix; };

    //! return covariance matrix
    EMatrix GetCovarianceMatrix() const { return *cov_matrix; };

    //! return the number of taxa
    int GetNtaxa() const { return Ntaxa; }

    //! return the exponentiated value of the multivariate brownian process for a given node and
    //! a given dimensions of the process
    double GetExpBrownianEntry(Tree::NodeIndex node, int dim) const {
        return exp(GetBrownianEntry(node, dim));
    }

    //! return the value of the multivariate brownian process for a given node and a given
    //! dimensions of the process
    double GetBrownianEntry(Tree::NodeIndex node, int dim) const {
        return node_multivariate->GetVal(node)(dim);
    }

    //! return number of dimensions of the multivariate brownian process
    int GetDimension() const { return dimensions; }

    std::string GetDimensionName(int dim) const {
        assert(taxon_traits != nullptr);
        return taxon_traits->GetHeader(dim);
    }


    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in Move), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void UpdateRoot() {
        for (int i = 0; i < dimensions; i++) {
            root_trait[i] = node_multivariate->GetVal(tree->root())(i);
        }
    }
    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() {
        chronogram->Update();
        chronogram->Scale();
        UpdateRoot();
    }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return PrecisionMatrixLogProb() + node_multivariate->GetLogProb(); }

    //! log prior of branch rate (brownian process) around of focal node
    double LocalNodeMultivariateLogPrior(Tree::NodeIndex node) const {
        return node_multivariate->GetLocalLogProb(node);
    }

    //! log prior of branch rate (brownian process) around of focal node
    double GetNodeLogProb(Tree::NodeIndex node) const {
        return node_multivariate->GetLogProb(node);
    }

    //! log prob of precision matrix
    double PrecisionMatrixLogProb() const {
        return prior_matrix->GetLogProb() + precision_matrix->GetLogProb(*prior_matrix);
    }

    double GetChronogramLength() const { return chronogram->GetSum(); }

    // Scatter (brownian process)
    void CollectScatterSuffStat() {
        scattersuffstat->Clear();
        scattersuffstat->AddSuffStat(*node_multivariate);
    }

    //! \brief complete MCMC move schedule
    double Move() {
        MoveParameters(30);
        return 1.0;
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveNodeTraits(0.5, 3, true);
            MoveNodeTraits(0.5, 3, false);
            MoveNodeTraits(0.05, 3, true);
            MoveNodeTraits(0.05, 3, false);

            CollectScatterSuffStat();
            SamplePrecisionMatrix();
            MovePriorMatrix(0.1, 3);
            MovePriorMatrix(0.01, 3);
        }
        precision_matrix->UpdateCovarianceMatrix(*cov_matrix);
        UpdateRoot();
    }

    void SamplePrecisionMatrix() {
        scattersuffstat->SamplePrecisionMatrix(*precision_matrix, *prior_matrix);
    };


    //! recompute the precision matrix and update the covariance matrix
    void RecomputePrecisionMatrix() {
        CollectScatterSuffStat();
        SamplePrecisionMatrix();
        precision_matrix->UpdateCovarianceMatrix(*cov_matrix);
    }

    //! MH moves on the invert wishart matrix (prior of the covariance matrix)
    void MovePriorMatrix(double tuning, int nrep) {
        if (uniq_kappa) {
            for (int rep = 0; rep < nrep; rep++) {
                double deltalogprob = -PrecisionMatrixLogProb();
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*prior_matrix) *= e;
                deltalogprob += PrecisionMatrixLogProb();
                deltalogprob += m;
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (!accepted) { (*prior_matrix) /= e; }
            }
        } else {
            for (int i = 0; i < dimensions; ++i) {
                Move::Scaling((*prior_matrix)(i), tuning, nrep,
                    &DatedNodeModel::PrecisionMatrixLogProb, &DatedNodeModel::NoUpdate, this);
            }
        }
    }

    //! MH moves on traits (brownian process)
    void MoveNodeTraits(double tuning, int nrep, bool leaves_to_root) {
        if (taxon_traits == nullptr or taxon_traits->GetDim() == 0) { return; }
        for (int rep = 0; rep < nrep; rep++) {
            for (int trait_dim = 0; trait_dim < taxon_traits->GetDim(); trait_dim++) {
                for (Tree::NodeIndex node :
                    leaves_to_root ? tree->leaves_root_to_iter() : tree->root_to_leaves_iter()) {
                    if (!tree->is_leaf(node) or
                        !taxon_traits->DataPresence(taxonmap->NodeToTaxon(node), trait_dim)) {
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
    }


    void ToStream(std::ostream &os) { os << *this; }
};

std::istream &operator>>(std::istream &is, std::unique_ptr<DatedNodeModel> &m) {
    std::string model_name, treefile, traitsfile;
    int prior_cov_df;
    bool uniq_kappa;
    is >> model_name;
    if (model_name != "DatedNode") {
        std::cerr << "Expected DatedNode for model name, got " << model_name << "\n";
        exit(1);
    }

    is >> treefile >> traitsfile >> prior_cov_df >> uniq_kappa;
    m = std::make_unique<DatedNodeModel>(treefile, traitsfile, prior_cov_df, uniq_kappa);
    Tracer tracer{*m};
    tracer.read_line(is);
    m->Update();
    return is;
}

std::ostream &operator<<(std::ostream &os, DatedNodeModel &m) {
    Tracer tracer{m};
    os << "DatedNode" << '\t';
    os << m.treefile << '\t';
    os << m.traitsfile << '\t';
    os << m.prior_cov_df << '\t';
    os << m.uniq_kappa << '\t';
    tracer.write_line(os);
    return os;
}
