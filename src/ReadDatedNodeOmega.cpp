#include <cmath>
#include <fstream>
#include "DatedNodeOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"
#include "tree/export.hpp"

using namespace std;
using namespace TCLAP;


class ReadNodeOmegaArgParse : public ReadArgParse {
  public:
    explicit ReadNodeOmegaArgParse(CmdLine &cmd) : ReadArgParse(cmd) {}

    SwitchArg newick{"t", "newick",
        "Computes the mean posterior node-specific entries of the multivariate Brownian process. "
        "Each entry of the multivariate Brownian process is written in a newick extended (.nhx) "
        "format file."
        "For each trait, results are written in {chain_name}.{trait}.nhx by default (optionally "
        "use the --output argument to specify a different output path).",
        cmd};
    SwitchArg newick_trees{"n", "newick_trees",
        "Export the node-specific entries of the multivariate Brownian process for each point of "
        "the MCMC. "
        "All entries of the multivariate Brownian process are written in a single file (.trees), "
        "containing as many lines as points in the MCMC. "
        "Each point of the MCMC (each line of the .trees file) is formatted as a newick extended "
        "tree (NHX). "
        "Results are written in {chain_name}.trees by default (optionally use the --output "
        "argument to specify a different output path).",
        cmd};
    SwitchArg same_space_as_input_traits{"l", "same_space_as_input_traits",
        "Export the node-specific entries of the multivariate Brownian process in the same space "
        "as in the input file provided by the option --traitsfile. "
        "By default the input values are assumed to be in log-space and hence the output are in "
        "the natural space, meaning the output are exponentiated value of the input."
        "This option makes sense if the trait that were provided in the input file were not "
        "intended to be transformed and you want them in the same space."
        "This option is only used with --newick and --newick_trees options.",
        cmd};
    SwitchArg cov{"c", "cov",
        "Computes the mean posterior covariance matrix, precision matrix and correlation matrix. "
        "Results are written in {chain_name}.cov by default (optionally use the --output argument "
        "to specify a different output path).",
        cmd};
};


int main(int argc, char *argv[]) {
    CmdLine cmd{"DatedMutSel", ' ', "0.1"};
    ReadNodeOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    unique_ptr<DatedNodeOmegaModel> model = nullptr;
    new ChainDriver(is);
    is >> model;
    ChainReader cr(*model, chain_name + ".chain");

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.GetPpred()) {
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model->PostPred("ppred_" + chain_name + "_" + to_string(i) + ".ali");
        }
        cerr << '\n';
    } else if (read_args.trace.getValue()) {
        string file_name = read_args.OutputFile(".trace");
        recompute_trace<DatedNodeOmegaModel>(*model, cr, file_name, every, size);
    } else if (read_args.cov.getValue()) {
        string file_name = read_args.OutputFile(".cov");
        ofstream os(file_name);
        os << "entries are in the following order:" << endl;

        for (int dim = 0; dim < model->GetDimension(); dim++) {
            os << model->GetDimensionName(dim) << endl;
        }
        EMatrix cov_matrix = EMatrix::Zero(model->GetDimension(), model->GetDimension());
        EMatrix posterior_prob = EMatrix::Zero(model->GetDimension(), model->GetDimension());
        EMatrix precision_matrix = EMatrix::Zero(model->GetDimension(), model->GetDimension());
        EMatrix partial_posterior_prob =
            EMatrix::Zero(model->GetDimension(), model->GetDimension());
        int count = 0;
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            EMatrix cov_matrix_chain = model->GetCovarianceMatrix();
            EMatrix precision_matrix_chain = model->GetPrecisionMatrix();
            bool is_positive = true;
            for (int i = 0; i < model->GetDimension(); i++) {
                if (cov_matrix_chain(i, i) <= 0) { is_positive = false; }
            }
            if (!is_positive) {
                model->RecomputePrecisionMatrix();
                cov_matrix_chain = model->GetCovarianceMatrix();
                precision_matrix_chain = model->GetPrecisionMatrix();
            }
            count++;
            cov_matrix += cov_matrix_chain;
            precision_matrix += precision_matrix_chain;
            for (int i = 0; i < model->GetDimension(); i++) {
                for (int j = 0; j < model->GetDimension(); j++) {
                    if (cov_matrix_chain.coeffRef(i, j) > 0) { posterior_prob.coeffRef(i, j) += 1; }
                    if (precision_matrix_chain.coeffRef(i, j) < 0) {
                        partial_posterior_prob.coeffRef(i, j) += 1;
                    }
                }
            }
        }
        cov_matrix /= count;
        posterior_prob /= count;
        precision_matrix /= count;
        partial_posterior_prob /= count;
        std::cerr << "Counted " << count << " steps out of " << size << endl;
        EMatrix cor_matrix = EMatrix::Zero(model->GetDimension(), model->GetDimension());
        EMatrix partial_cor_matrix = EMatrix::Zero(model->GetDimension(), model->GetDimension());
        for (int i = 0; i < model->GetDimension(); i++) {
            for (int j = 0; j < model->GetDimension(); j++) {
                cor_matrix.coeffRef(i, j) =
                    cov_matrix.coeffRef(i, j) /
                    sqrt(cov_matrix.coeffRef(i, i) * cov_matrix.coeffRef(j, j));
                partial_cor_matrix.coeffRef(i, j) =
                    -precision_matrix.coeffRef(i, j) /
                    sqrt(precision_matrix.coeffRef(i, i) * precision_matrix.coeffRef(j, j));
            }
        }
        export_matrix(os, model->GetDimension(), cov_matrix, "covariances");
        export_matrix(os, model->GetDimension(), cor_matrix, "correlation coefficients");
        export_matrix(os, model->GetDimension(), posterior_prob,
            "posterior probabilities of a positive coefficient", false);
        export_matrix(os, model->GetDimension(), precision_matrix, "precisions");
        export_matrix(
            os, model->GetDimension(), partial_cor_matrix, "partial correlation coefficients");
        export_matrix(os, model->GetDimension(), partial_posterior_prob,
            "posterior probabilities of a positive partial coefficient", false);

        cerr << endl << "matrices in " << file_name << "." << endl;
    } else if (read_args.newick.getValue()) {
        vector<vector<vector<double>>> dim_node_traces(model->GetDimension());
        vector<vector<double>> branch_times(model->GetTree().nb_nodes());
        vector<vector<double>> branch_length(model->GetTree().nb_nodes());
        vector<vector<double>> branch_omega(model->GetTree().nb_nodes());

        for (int dim{0}; dim < model->GetDimension(); dim++) {
            dim_node_traces[dim].resize(model->GetTree().nb_nodes());
        }

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);

            model->UpdateBranches(true);
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
                node++) {
                if (!model->GetTree().is_root(node)) {
                    double branch_time = model->GetBranchTime(node);
                    assert(branch_time >= 0);
                    branch_times[node].push_back(branch_time);
                    branch_length[node].push_back(model->GetBranchLength(node));
                    branch_omega[node].push_back(model->GetBranchOmega(node));
                }
                for (int dim{0}; dim < model->GetDimension(); dim++) {
                    if (read_args.same_space_as_input_traits.getValue()) {
                        dim_node_traces[dim][node].push_back(model->GetBrownianEntry(node, dim));
                    } else {
                        dim_node_traces[dim][node].push_back(model->GetExpBrownianEntry(node, dim));
                    }
                }
            }
        }
        cerr << '\n';

        ExportTree base_export_tree(model->GetTree());
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
            node++) {
            if (!model->GetTree().is_root(node)) {
                base_export_tree.set_tag(node, "length", to_string(mean(branch_times[node])));
            }
        }

        export_tree(base_export_tree, "BranchLength", read_args.OutputFile(), branch_length);
        export_tree(base_export_tree, "BranchTime", read_args.OutputFile(), branch_times);
        export_tree(base_export_tree, "BranchOmega", read_args.OutputFile(), branch_omega);
        for (int dim{0}; dim < model->GetDimension(); dim++) {
            export_tree(base_export_tree, model->GetDimensionName(dim), read_args.OutputFile(),
                dim_node_traces[dim]);
        }
    } else if (read_args.newick_trees.getValue()) {
        string nhxname = read_args.OutputFile(".trees");
        std::ofstream nhx(nhxname);
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);

            model->UpdateBranches(true);
            ExportTree export_tree(model->GetTree());
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
                node++) {
                if (!model->GetTree().is_root(node)) {
                    export_tree.set_tag(node, "length", to_string(model->GetBranchTime(node)));
                    export_tree.set_tag(
                        node, "BranchLength", to_string(model->GetBranchLength(node)));
                    export_tree.set_tag(
                        node, "BranchOmega", to_string(model->GetBranchOmega(node)));
                    export_tree.set_tag(node, "BranchMutationRatePerTime",
                        to_string(model->GetBranchMutRate(node)));
                }
                for (int dim{0}; dim < model->GetDimension(); dim++) {
                    if (read_args.same_space_as_input_traits.getValue()) {
                        export_tree.set_tag(node, model->GetDimensionName(dim),
                            to_string(model->GetBrownianEntry(node, dim)));
                    } else {
                        export_tree.set_tag(node, model->GetDimensionName(dim),
                            to_string(model->GetExpBrownianEntry(node, dim)));
                    }
                }
            }
            nhx << export_tree.as_string() << std::endl;
        }
        cerr << '\n';
        nhx.close();
        std::cerr << "Trees in " << nhxname << "\n";
    } else {
        stats_posterior<DatedNodeOmegaModel>(*model, cr, every, size);
    }
}