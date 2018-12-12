#include <cmath>
#include <fstream>
#include "DatedBranchOmegaModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "DatedBranchOmega", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DatedBranchOmegaModel> model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        is >> model;
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = unique_ptr<DatedBranchOmegaModel>(
            new DatedBranchOmegaModel(args.alignment.getValue(), args.treefile.getValue()));
    }

    model->Update();
    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();

    delete chain_driver;
}
