#include <cmath>
#include <fstream>
#include "DatedNodeMutSelModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"
#include "tree/export.hpp"

using namespace std;
using namespace TCLAP;

std::string double2str(double val) {
    std::ostringstream so;
    so << std::scientific << val;
    return so.str();
}

class ReadNodeMutSelArgParse : public ReadArgParse {
  public:
    explicit ReadNodeMutSelArgParse(CmdLine &cmd) : ReadArgParse(cmd) {}

    TCLAP::ValueArg<string> profiles{"o", "profiles",
        "Output profiles name if desired (otherwise given by {chain_name}.siteprofiles)", false, "",
        "string", cmd};

    SwitchArg ss{
        "s", "ss", "Computes the mean posterior site-specific state equilibrium frequencies", cmd};

    SwitchArg newick{"t", "newick",
        "Computes the mean posterior node-specific entries of the multivariate Brownian process",
        cmd};

    string GetProfilesName() {
        if (profiles.getValue().empty()) {
            return GetChainName() + ".siteprofiles";
        } else {
            return profiles.getValue();
        }
    }
};

void export_tree(ExportTree export_tree, string const &name, string const &path,
    vector<vector<double>> &values) {
    for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(export_tree.GetTree().nb_nodes());
         node++) {
        if (!values[node].empty()) {
            sort(values[node].begin(), values[node].end());
            auto up = static_cast<size_t>(0.95 * values[node].size());
            if (up >= values.size()) { up = values[node].size() - 1; }
            auto down = static_cast<size_t>(0.05 * values[node].size());
            export_tree.set_tag(node, name + "_min", double2str(values[node].at(down)));
            export_tree.set_tag(node, name + "_max", double2str(values[node].at(up)));
            export_tree.set_tag(node, name, double2str(mean(values[node])));
        }
    }
    string nhxname = path + "." + name + ".nhx";
    ofstream nhx(nhxname);
    nhx << export_tree.as_string() << endl;
    nhx.close();
    cerr << "Tree for " << name << " in " << nhxname << "\n";
}

int main(int argc, char *argv[]) {
    CmdLine cmd{"DatedMutSel", ' ', "0.1"};
    ReadNodeMutSelArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    ChainDriver *fake_read = nullptr;
    unique_ptr<DatedNodeMutSelModel> model = nullptr;
    fake_read = new ChainDriver(is);
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
        recompute_trace<DatedNodeMutSelModel>(*model, cr, chain_name, every, size);
    } else if (read_args.ss.getValue()) {
        vector<vector<double>> sitestat(model->GetNsite(), {0});

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model->GetNsite(); i++) {
                vector<double> const &profile = model->GetProfile(i);
                if (sitestat[i].size() != profile.size()) { sitestat[i].resize(profile.size(), 0); }
                for (unsigned k{0}; k < profile.size(); k++) { sitestat[i][k] += profile[k]; }
            }
        }
        cerr << '\n';

        ofstream os(read_args.GetProfilesName().c_str());
        os << model->GetNsite() << '\n';
        for (int i = 0; i < model->GetNsite(); i++) {
            os << i + 1;
            for (auto &aa : sitestat[i]) {
                aa /= size;
                os << '\t' << aa;
            }
            os << '\n';
        }
        cerr << "mean site-specific profiles in " << read_args.GetProfilesName() << "\n";
        cerr << '\n';
    } else if (read_args.newick.getValue()) {
        vector<vector<vector<double>>> dim_node_traces(model->GetDimension());
        vector<vector<double>> branch_times(model->GetTree().nb_nodes());
        vector<vector<double>> log10_branch_length(model->GetTree().nb_nodes());
        vector<vector<double>> log10_leaves_theta(model->GetTree().nb_nodes());
        vector<vector<double>> contrast_pop_size(model->GetTree().nb_nodes());

        for (int dim{0}; dim < model->GetDimension(); dim++) {
            dim_node_traces[dim].resize(model->GetTree().nb_nodes());
        }

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);

            model->UpdateBranches();
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
                 node++) {
                if (!model->GetTree().is_root(node)) {
                    double branch_time = model->GetBranchTime(node);
                    assert(branch_time >= 0);
                    assert(branch_time <= 1);
                    branch_times[node].push_back(branch_time);
                    log10_branch_length[node].push_back(log10(model->GetBranchLength(node)));
                    contrast_pop_size[node].push_back(model->GetContrast(node, dim_pop_size));
                }
                if (model->PolymorphismAware() and model->GetTree().is_leaf(node)) {
                    log10_leaves_theta[node].push_back(log10(model->GetTheta(node)));
                }
                for (int dim{0}; dim < model->GetDimension(); dim++) {
                    dim_node_traces[dim][node].push_back(model->GetBrownianEntry(node, dim));
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

        export_tree(base_export_tree, "ContrastPopulationSize", read_args.GetChainName(),
            contrast_pop_size);
        export_tree(
            base_export_tree, "Log10BranchLength", read_args.GetChainName(), log10_branch_length);
        export_tree(base_export_tree, "BranchTime", read_args.GetChainName(), branch_times);
        if (model->PolymorphismAware()) {
            export_tree(
                base_export_tree, "Log10Theta", read_args.GetChainName(), log10_leaves_theta);
        }
        for (int dim{0}; dim < model->GetDimension(); dim++) {
            export_tree(base_export_tree, model->GetDimensionName(dim), read_args.GetChainName(),
                dim_node_traces[dim]);
        }
    } else {
        stats_posterior<DatedNodeMutSelModel>(*model, cr, every, size);
    }
}