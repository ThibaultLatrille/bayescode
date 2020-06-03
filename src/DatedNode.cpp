#include <cmath>
#include <fstream>
#include "DatedNodeModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class DatedNodeOmegaArgParse : public BaseArgParse {
  public:
    explicit DatedNodeOmegaArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    ValueArg<std::string> traitsfile{
        "", "traitsfile", "Traits file for taxon at the leaves", false, "Null", "string", cmd};
    ValueArg<std::string> fossils{
        "", "fossils", "Fossils data (to clamp the node ages)", false, "Null", "string", cmd};
    ValueArg<int> prior_cov_df{"", "df", "Invert Wishart degree of freedom", false, 0, "int", cmd};
    SwitchArg uniq_kappa{"", "uniq_kappa",
        "Unique kappa for the invert Wishart matrix prior (otherwise 1 for each dimension)", cmd,
        false};
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "DatedNodeOmega", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DatedNodeModel> model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        is >> model;
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        TreeAppArgParse inference_args(cmd);
        DatedNodeOmegaArgParse args(cmd);
        cmd.parse();
        chain_driver = new ChainDriver(
            cmd.chain_name(), inference_args.every.getValue(), inference_args.until.getValue());
        model = std::make_unique<DatedNodeModel>(inference_args.treefile.getValue(),
            args.traitsfile.getValue(), args.fossils.getValue(), args.prior_cov_df.getValue(),
            args.uniq_kappa.getValue());
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

    delete chain_driver;
}
