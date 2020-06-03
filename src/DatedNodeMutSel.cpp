#include <cmath>
#include <fstream>
#include "DatedNodeMutSelModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class DatedNodeMutselArgParse : public BaseArgParse {
  public:
    explicit DatedNodeMutselArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    ValueArg<int> ncat{
        "", "ncat", "truncation of the first-level stick-breaking process", false, 100, "int", cmd};
    ValueArg<int> basencat{"", "basencat", "truncation of the second-level stick-breaking process",
        false, 1, "int", cmd};
    SwitchArg condition_aware{"b", "condition_aware",
        "One Ne per condition, if the tree doesn't have condition, then one Ne per branch", cmd,
        false};
    ValueArg<std::string> traitsfile{
        "", "traitsfile", "Traits file for taxon at the leaves", false, "Null", "string", cmd};
    ValueArg<std::string> profiles{
        "c", "profiles", "Preferences profiles (to clamp)", false, "Null", "string", cmd};
    SwitchArg move_root_pop_size{"", "move_root_pop_size",
        "Move Ne at the root (for the equilibrium frequencies)", cmd, false};
    SwitchArg clamp_pop_sizes{
        "", "clamp_pop_sizes", "Clamp the branch population size", cmd, false};
    SwitchArg clamp_nuc_matrix{"", "clamp_nuc_matrix", "Clamp the nucleotide matrix", cmd, false};
    SwitchArg clamp_corr_matrix{
        "", "clamp_corr_matrix", "Clamp the correlation matrix", cmd, false};
    SwitchArg polymorphism_aware{"p", "polymorphism_aware", "Use polymorphic data", cmd, false};
    ValueArg<unsigned> precision{
        "", "precision", "The precision of PRF computation", false, 6, "unsigned", cmd};
    SwitchArg arithmetic{
        "d", "arithmetic", "Use arithmetic mean instead of arithmetic", cmd, false};
    ValueArg<std::string> fossils{
        "", "fossils", "Fossils data (to clamp the node ages)", false, "Null", "string", cmd};
    ValueArg<int> prior_cov_df{"", "df", "Invert Wishart degree of freedom", false, 0, "int", cmd};
    SwitchArg uniq_kappa{"", "uniq_kappa",
        "Unique kappa for the invert Wishart matrix prior (otherwise 1 for each dimension)", cmd,
        false};
    void check() {
        if (condition_aware.getValue()) {
            cerr << "The switch parameter ([-b] or [--condition_aware]) is not yet implemented."
                 << endl;
        }
        if (profiles.getValue() != "Null") {
            cout << "Preferences are clamped (option [--profiles <string>]), thus options [--ncat] "
                    "and [-basencat] are not used"
                 << endl;
        }
    }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "DatedMutSel", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DatedNodeMutSelModel> model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        is >> model;
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse inference_args(cmd);
        DatedNodeMutselArgParse args(cmd);
        cmd.parse();
        args.check();
        chain_driver = new ChainDriver(
            cmd.chain_name(), inference_args.every.getValue(), inference_args.until.getValue());
        model = std::make_unique<DatedNodeMutSelModel>(inference_args.alignment.getValue(),
            inference_args.treefile.getValue(), args.traitsfile.getValue(),
            args.profiles.getValue(), args.ncat.getValue(), args.basencat.getValue(),
            args.condition_aware.getValue(), args.polymorphism_aware.getValue(),
            args.precision.getValue(), args.arithmetic.getValue(),
            args.move_root_pop_size.getValue(), args.clamp_pop_sizes.getValue(),
            args.clamp_nuc_matrix.getValue(), args.clamp_corr_matrix.getValue(),
            args.fossils.getValue(), args.prior_cov_df.getValue(), args.uniq_kappa.getValue());
        model->Update();
    }
    model->ResampleSub(1.0);
    model->MoveParameters(10);
    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();
}
