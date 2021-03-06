#include <cmath>
#include <fstream>
#include "DatedNodeModel.hpp"
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

class ReadNodeOmegaArgParse : public ReadArgParse {
  public:
    explicit ReadNodeOmegaArgParse(CmdLine &cmd) : ReadArgParse(cmd) {}

    SwitchArg newick{"t", "newick",
        "Computes the mean posterior node-specific entries of the multivariate Brownian process",
        cmd};
};

void export_tree(ExportTree export_tree, string name, string const &path,
    vector<vector<double>> &values, bool log_value = true) {
    if (log_value) { name = "Log" + name; }
    for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(export_tree.GetTree().nb_nodes());
         node++) {
        auto value_array = values[node];
        if (!log_value) {
            for (auto &v : value_array) { v = exp(v); }
        }
        if (!value_array.empty()) {
            sort(value_array.begin(), value_array.end());
            auto up = static_cast<size_t>(0.95 * value_array.size());
            if (up >= values.size()) { up = value_array.size() - 1; }
            auto down = static_cast<size_t>(0.05 * value_array.size());
            export_tree.set_tag(node, name + "_min", double2str(value_array.at(down)));
            export_tree.set_tag(node, name + "_max", double2str(value_array.at(up)));
            export_tree.set_tag(node, name, double2str(mean(value_array)));
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
    ReadNodeOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    ChainDriver *fake_read = nullptr;
    unique_ptr<DatedNodeModel> model = nullptr;
    fake_read = new ChainDriver(is);
    is >> model;
    ChainReader cr(*model, chain_name + ".chain");

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.trace.getValue()) {
        recompute_trace<DatedNodeModel>(*model, cr, chain_name, every, size);
    } else if (read_args.newick.getValue()) {
        vector<vector<vector<double>>> dim_node_traces(model->GetDimension());
        vector<vector<double>> branch_times(model->GetTree().nb_nodes());

        for (int dim{0}; dim < model->GetDimension(); dim++) {
            dim_node_traces[dim].resize(model->GetTree().nb_nodes());
        }

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);

            model->Update();
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
                 node++) {
                if (!model->GetTree().is_root(node)) {
                    double branch_time = model->GetBranchTime(node);
                    assert(branch_time >= 0);
                    assert(branch_time <= 1);
                    branch_times[node].push_back(branch_time);
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

        export_tree(base_export_tree, "BranchTime", read_args.GetChainName(), branch_times, false);
        for (int dim{0}; dim < model->GetDimension(); dim++) {
            export_tree(base_export_tree, model->GetDimensionName(dim), read_args.GetChainName(),
                dim_node_traces[dim], true);
            export_tree(base_export_tree, model->GetDimensionName(dim), read_args.GetChainName(),
                dim_node_traces[dim], false);
        }
    } else {
        stats_posterior<DatedNodeModel>(*model, cr, every, size);
    }
}