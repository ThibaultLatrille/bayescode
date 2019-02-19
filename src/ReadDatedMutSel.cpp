#include <cmath>
#include <fstream>
#include "DatedMutSelModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;

class ReadDatedMutSelArgParse : public ReadArgParse {
  public:
    explicit ReadDatedMutSelArgParse(CmdLine &cmd) : ReadArgParse(cmd) {}

    TCLAP::ValueArg<std::string> profiles{"o", "profiles",
        "Output profiles name if desired (otherwise given by {chain_name}.siteprofiles)", false, "", "string", cmd};

    SwitchArg ss{
        "s", "ss", "Computes the mean posterior site-specific state equilibrium frequencies", cmd};

    std::string GetProfilesName() {
        if (profiles.getValue().empty()) {
            return GetChainName() + ".siteprofiles";
        } else {
            return profiles.getValue();
        }
    }
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"DatedMutSel", ' ', "0.1"};
    ReadDatedMutSelArgParse read_args(cmd);
    cmd.parse(argc, argv);

    std::string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    std::ifstream is{chain_name + ".param"};
    ChainDriver *fake_read = nullptr;
    unique_ptr<DatedMutSelModel> model = nullptr;
    fake_read = new ChainDriver(is);
    is >> model;
    ChainReader cr(*model, chain_name + ".chain");

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.GetPpred()) {
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model->PostPred(
                "ppred_" + chain_name + "_" + std::to_string(i) + ".ali");
        }
        cerr << '\n';
    } else if (read_args.ss.getValue()) {
        std::vector<std::vector<double>> sitestat(model->GetNsite(), {0});

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model->GetNsite(); i++) {
                std::vector<double> const &profile = model->GetProfile(i);
                if (sitestat[i].size() != profile.size()) {
                    sitestat[i].resize(profile.size(), 0);
                };
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
    } else {
        stats_posterior<DatedMutSelModel>(*model, cr, every, size);
    }
}