#include <cmath>
#include <fstream>
#include "DatedMutSelModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class DatedMutselArgParse : public BaseArgParse {
  public:
    DatedMutselArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    ValueArg<int> ncat{
        "", "ncat", "truncation of the first-level stick-breaking process", false, 100, "int", cmd};
    ValueArg<int> basencat{"", "basencat", "truncation of the second-level stick-breaking process",
        false, 1, "int", cmd};
    SwitchArg condition_aware{"b", "condition_aware",
        "One Ne per condition, if the tree doesn't have condition, then one Ne per branch", cmd,
        false};
    ValueArg<unsigned> precision{"", "precision", "The precision of PRF computation", false, 6, "unsigned", cmd};
    ValueArg<std::string> profiles{"c", "profiles", "Preferences profiles (to clamp)", false, "", "string", cmd};
    SwitchArg polymorphism_aware{"p", "polymorphism_aware", "Use polymorphic data", cmd, false};
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "DatedMutSel", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DatedMutSelModel> model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        is >> model;
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse args(cmd);
        DatedMutselArgParse datedmutsel_args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = unique_ptr<DatedMutSelModel>(new DatedMutSelModel(args.alignment.getValue(),
            args.treefile.getValue(), datedmutsel_args.profiles.getValue(),
            datedmutsel_args.ncat.getValue(), datedmutsel_args.basencat.getValue(),
            datedmutsel_args.condition_aware.getValue(), datedmutsel_args.polymorphism_aware.getValue(), datedmutsel_args.precision.getValue()));
        model->Update();
    }

    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();
}
