#pragma once

#include "Chronogram.hpp"
#include "Move.hpp"
#include "MultivariateProcess.hpp"
#include "ScatterSuffStat.hpp"
#include "TaxonTraits.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"

/**
 * \brief Traits are modeled as a browian process along the tree (like CoEvol)
 *
 */


class DatedNodeModel : public ChainComponent {
    // tree and data
    std::string treefile, traitsfile{"Null"}, fossilsfile{"Null"};
    std::unique_ptr<const Tree> tree;
    TaxonSet *taxonset{nullptr};
    TaxonMap *taxonmap{nullptr};
    TaxonTraits *taxon_traits{nullptr};

    int Ntaxa;
    int Nbranch;

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
    ScatterSuffStat *scattersuffstat;

  public:
    friend std::ostream &operator<<(std::ostream &os, DatedNodeModel &m);

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    DatedNodeModel(std::string intreefile, std::string intraitsfile, std::string infossilsfile,
        int inprior_cov_df, bool inuniq_kappa)
        : treefile(std::move(intreefile)),
          traitsfile(std::move(intraitsfile)),
          fossilsfile(std::move(infossilsfile)),
          prior_cov_df{inprior_cov_df},
          uniq_kappa{inuniq_kappa} {
        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);
        Nbranch = tree->nb_branches();

        // Node ages
        nodeages = new NodeAges(*tree, fossilsfile);
        chronogram = new Chronogram(*nodeages);

        std::vector<std::string> leave_names;
        for (auto const &n : tree->root_to_leaves_iter()) {
            if (tree->is_leaf(n)) { leave_names.push_back(tree->node_name(n)); }
        }
        taxonset = new TaxonSet(leave_names);
        taxonmap = new TaxonMap(tree.get(), taxonset);
        Ntaxa = taxonmap->GetNtaxa();
        if (traitsfile != "Null") {
            taxon_traits = new TaxonTraits(traitsfile, *taxonset, false, 0);
        }

        // Precision matrix and prior
        dimensions = taxon_traits->GetDim();
        prior_matrix = new PriorCovariance(dimensions, prior_cov_df, uniq_kappa);
        precision_matrix = new PrecisionMatrix(*prior_matrix);

        node_multivariate = new NodeMultivariateProcess(*chronogram, *precision_matrix, dimensions);
        if (taxon_traits != nullptr) { node_multivariate->ClampLeaves(*taxon_traits, *taxonmap); }
        // Suff Stats
        scattersuffstat = new ScatterSuffStat(*tree);
    }

    virtual ~DatedNodeModel() = default;

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        model_stat(info, "lnprob", [&]() { return DatedNodeModel::GetLogProb(); });
        model_node(info, "nodeages", *nodeages);
        model_node(info, "node_multivariate", *node_multivariate);
        model_node(info, "prior_cov_matrix", *prior_matrix);
        model_node(info, "precision_matrix", *precision_matrix);

        for (int i = 0; i < dimensions; i++) {
            model_stat(info, "PriorSigma_" + std::to_string(i), prior_matrix->coeffRef(i));
            for (int j = 0; j <= i; j++) {
                model_stat(info, "Precision_" + std::to_string(i) + "_" + std::to_string(j),
                    precision_matrix->coeffRef(i, j));
            }
        }
        // Descriptive statistics - for each branch of the tree
        for (Tree::BranchIndex branch = 0; branch < tree->nb_branches(); branch++) {
            string b_name = tree->node_name(tree->node_index(branch));
            model_stat(info, "*BranchTime_" + b_name, (*chronogram)[branch]);
        }
    }

    //! return tree
    const Tree &GetTree() const { return *tree; }

    //! return time of branch
    double GetBranchTime(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return chronogram->GetVal(tree->branch_index(node));
    };

    //! return the value of the multivariate brownian process for a given node and a given
    //! dimensions of the process
    double GetBrownianEntry(Tree::NodeIndex node, int dim) const {
        return node_multivariate->GetVal(node)(dim);
    }

    //! return number of dimensions of the multivariate brownian process
    int GetDimension() const { return dimensions; }

    std::string GetDimensionName(int dim) const {
        assert(taxon_traits != nullptr);
        return "Traits" + taxon_traits->GetHeader(dim);
    }


    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in Move), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() { chronogram->Update(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return PrecisionMatrixLogProb() + node_multivariate->GetLogProb(); }

    //! log prior of branch rate (brownian process) around of focal node
    double LocalNodeMultivariateLogPrior(Tree::NodeIndex node) const {
        return node_multivariate->GetLocalLogProb(node);
    }

    //! log prob of precision matrix
    double PrecisionMatrixLogProb() const {
        return prior_matrix->GetLogProb() + precision_matrix->GetLogProb(*prior_matrix);
    }

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
    }

    void SamplePrecisionMatrix() {
        scattersuffstat->SamplePrecisionMatrix(*precision_matrix, *prior_matrix);
    };

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
    std::string model_name, treefile, traitsfile, fossilsfile;
    int prior_cov_df;
    bool uniq_kappa;
    is >> model_name;
    if (model_name != "DatedNode") {
        std::cerr << "Expected DatedNode for model name, got " << model_name << "\n";
        exit(1);
    }

    is >> treefile >> traitsfile >> fossilsfile >> prior_cov_df >> uniq_kappa;
    m = std::make_unique<DatedNodeModel>(
        treefile, traitsfile, fossilsfile, prior_cov_df, uniq_kappa);
    Tracer tracer{*m};
    tracer.read_line(is);
    m->Update();
    return is;
}

std::ostream &operator<<(std::ostream &os, DatedNodeModel &m) {
    Tracer tracer{m};
    os << "DatedNode" << '\t';
    os << m.treefile << '\t';
    assert(!m.traitsfile.empty());
    os << m.traitsfile << '\t';
    assert(!m.fossilsfile.empty());
    os << m.fossilsfile << '\t';
    os << m.prior_cov_df << '\t';
    os << m.uniq_kappa << '\t';
    tracer.write_line(os);
    return os;
}
