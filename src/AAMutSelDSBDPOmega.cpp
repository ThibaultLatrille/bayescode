#include <cmath>
#include <fstream>
#include "AAMutSelDSBDPOmegaModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class AAMutselArgParse : public BaseArgParse {
  public:
    AAMutselArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    ValueArg<int> ncat{
        "", "ncat", "truncation of the first-level stick-breaking process", false, 100, "int", cmd};
    ValueArg<int> basencat{"", "basencat", "truncation of the second-level stick-breaking process",
        false, 1, "int", cmd};
    SwitchArg freeomega{"", "freeomega",
        "omega is allowed to vary with shrinkage (otherwise set to 1)", cmd, false};
    SwitchArg mixomega{"", "mixomega",
        "the prior of omega is a mix of (1-pi) at 1 and pi at 1+d, "
        "with d~Gamma(dposomhypermean,dposomhyperinvshape)",
        cmd, false};
    ValueArg<double> dposompi{"", "dposompi",
        "proportion of positively selected sites (if -mixomega is switched on)", false, 0.1,
        "double", cmd};
    ValueArg<double> dposomhypermean{"", "dposomhypermean",
        "inverse shape of the gamma distribution "
        "(if --mixomega is switched on)",
        false, 1.0, "double", cmd};
    ValueArg<double> dposomhyperinvshape{"", "dposomhyperinvshape",
        "inverse shape of the gamma distribution "
        "(if --mixomega is switched on)",
        false, 0.5, "double", cmd};
    SwitchArg polymorphism_aware{"p", "polymorphism_aware", "Use polymorphic data", cmd, false};
    ValueArg<unsigned> precision{
        "", "precision", "The precision of PRF computation", false, 6, "unsigned", cmd};

    //! - omegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    int omegamode() {
        if (freeomega.getValue()) {
            return 1;
        } else {
            return 3;
        }
    }

    int omegaprior() {
        if (mixomega.getValue()) {
            return 1;
        } else {
            return 0;
        }
    }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "AAMutSelDSBDPOmega", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    AAMutSelDSBDPOmegaModel *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        model = new AAMutSelDSBDPOmegaModel(is);
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse args(cmd);
        AAMutselArgParse aamutsel_args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = new AAMutSelDSBDPOmegaModel(args.alignment.getValue(), args.treefile.getValue(),
            aamutsel_args.omegamode(), aamutsel_args.omegaprior(),
            aamutsel_args.dposompi.getValue(), aamutsel_args.dposomhypermean.getValue(),
            aamutsel_args.dposomhyperinvshape.getValue(), aamutsel_args.ncat.getValue(),
            aamutsel_args.basencat.getValue(), aamutsel_args.polymorphism_aware.getValue(),
            aamutsel_args.precision.getValue());
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
