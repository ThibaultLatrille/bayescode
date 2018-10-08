#include <cmath>
#include <fstream>
#include "AAMutSelDSBDPOmegaModel.hpp"
#include "Chain.hpp"
#include "InferenceAppArgParse.hpp"
#include "components/ChainDriver.hpp"

using namespace std;

/**
 * \brief Chain object for running an MCMC under AAMutSelDSBDPOmegaModel
 */

class AAMutSelDSBDPOmegaChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int omegamode, omegaprior;
    double dposompi, dposomhypermean, dposomhyperinvshape;
    int Ncat, baseNcat;

  public:
    AAMutSelDSBDPOmegaModel *GetModel() { return static_cast<AAMutSelDSBDPOmegaModel *>(model); }

    string GetModelType() override { return modeltype; }

    AAMutSelDSBDPOmegaChain(string indatafile, string intreefile, int inomegamode, int inomegaprior,
        double indposompi, double indposomhypermean, double indposomhyperinvshape, int inNcat,
        int inbaseNcat, int inevery, int inuntil, string inname, int force)
        : modeltype("AAMUTSELDSBDPOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          omegamode(inomegamode),
          omegaprior(inomegaprior),
          dposompi(indposompi),
          dposomhypermean(indposomhypermean),
          dposomhyperinvshape(indposomhyperinvshape),
          Ncat(inNcat),
          baseNcat(inbaseNcat) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    AAMutSelDSBDPOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        cerr << "new model\n";
        model =
            new AAMutSelDSBDPOmegaModel(datafile, treefile, omegamode, omegaprior, Ncat, baseNcat);
        if (omegaprior == 1) {
            GetModel()->SetDPosOmHyperParameters(dposompi, dposomhypermean, dposomhyperinvshape);
        }
        cerr << "allocate\n";
        GetModel()->Allocate();
        cerr << "update\n";
        GetModel()->Update();
        cerr << "-- Reset" << endl;
        Reset(force);
        cerr << "-- initial ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        is >> omegamode >> omegaprior >> dposompi >> dposomhypermean >> dposomhyperinvshape >>
            Ncat >> baseNcat;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "AAMUTSELDSBDPOMEGA") {
            model = new AAMutSelDSBDPOmegaModel(
                datafile, treefile, omegamode, omegaprior, Ncat, baseNcat);
            if (omegaprior == 1) {
                GetModel()->SetDPosOmHyperParameters(
                    dposompi, dposomhypermean, dposomhyperinvshape);
            }
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        GetModel()->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << omegamode << '\t' << omegaprior << '\t' << dposompi << '\t' << dposomhypermean
                 << '\t' << dposomhyperinvshape << '\n';
        param_os << Ncat << '\t' << baseNcat << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

class AAMutselArgParse : public BaseArgParse {
  public:
    AAMutselArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}
    ValueArg<int> ncat{
        "", "ncat", "Maximum number of Dirichlet amino acid profiles", false, 100, "int", cmd};
    ValueArg<int> basencat{"", "basencat", "TODO", false, 1, "int", cmd};
    // SwitchArg fixomega{"", "fixomega", "TODO", cmd};
    SwitchArg freeomega{"", "freeomega", "TODO", cmd};
    SwitchArg gamomega{"", "gamomega", "TODO", cmd};
    SwitchArg mixomega{"", "mixomega", "TODO", cmd};
    SwitchArg omegaprior{"", "omegaprior", "TODO", cmd};
    ValueArg<double> dposompi{"", "dposompi", "TODO", false, 0.1, "double", cmd};
    ValueArg<double> dposomhypermean{"", "dposomhypermean", "TODO", false, 1.0, "double", cmd};
    ValueArg<double> dposomhyperinvshape{
        "", "dposomhyperinvshape", "TODO", false, 0.5, "double", cmd};

    int omegamode() {
        if (freeomega.getValue()) {
            return 1;
        } else {
            return 3;
        }
    }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "aamutsel", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    AAMutSelDSBDPOmegaModel *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        // model = new AAMutSelDSBDPOmegaModel(is); // TODO
    } else {
        InferenceAppArgParse args(cmd);
        AAMutselArgParse aamutsel_args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = new AAMutSelDSBDPOmegaModel(args.alignment.getValue(), args.treefile.getValue(),
            aamutsel_args.omegamode(), aamutsel_args.omegaprior.getValue(),
            aamutsel_args.ncat.getValue(), aamutsel_args.basencat.getValue());
    }


    string name = "";
    AAMutSelDSBDPOmegaChain *chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new AAMutSelDSBDPOmegaChain(name);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int Ncat = 100;
        int baseNcat = 1;
        int omegamode = 3;
        int omegaprior = 0;
        double dposompi = 0.1;
        double dposomhypermean = 1.0;
        double dposomhyperinvshape = 0.5;
        name = "";
        int force = 1;
        int every = 1;
        int until = -1;

        try {
            if (argc == 1) { throw(0); }

            int i = 1;
            while (i < argc) {
                string s = argv[i];

                if (s == "-d") {
                    i++;
                    datafile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-ncat") {
                    i++;
                    Ncat = atoi(argv[i]);
                } else if (s == "-basencat") {
                    i++;
                    baseNcat = atoi(argv[i]);
                } else if (s == "-fixomega") {
                    omegamode = 3;
                } else if (s == "-freeomega") {
                    omegamode = 1;
                } else if (s == "-gamomega") {
                    omegaprior = 0;
                } else if (s == "-mixomega") {
                    omegaprior = 1;
                    i++;
                    dposompi = atof(argv[i]);
                    i++;
                    dposomhypermean = atof(argv[i]);
                    i++;
                    dposomhyperinvshape = atof(argv[i]);
                } else if ((s == "-x") || (s == "-extract")) {
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                } else {
                    if (i != (argc - 1)) { throw(0); }
                    name = argv[i];
                }
                i++;
            }
            if ((datafile == "") || (treefile == "") || (name == "")) { throw(0); }
        } catch (...) {
            cerr << "aamutseldp -d <alignment> -t <tree> -ncat <ncat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new AAMutSelDSBDPOmegaChain(datafile, treefile, omegamode, omegaprior, dposompi,
            dposomhypermean, dposomhyperinvshape, Ncat, baseNcat, every, until, name, force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize()
         << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
