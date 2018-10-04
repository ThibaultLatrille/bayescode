#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneSingleOmegaModel.hpp"

using namespace std;

// MPI_Datatype Propagate_arg;

/**
 * \brief Chain object for running an MCMC under MultiGeneSingleOmegaModel
 */

class MultiGeneSingleOmegaChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int blmode, nucmode, omegamode;
    double omegahypermean, omegahyperinvshape;
    MultiGeneSingleOmegaModel* model;
    int myid;
    int nprocs;

    //! saving frequency (i.e. number of move cycles performed between each point
    //! saved to file)
    int every;
    //! intended final size of the chain (until==-1 means no a priori specified
    //! upper limit)
    int until;
    //! current size (number of points saved to file)
    int size;
    //! base name for all files corresponding to that chain
    string name;
    //! flag: if 1, then complete state is saved at each interation in .chain file
    int saveall;

public:
    MultiGeneSingleOmegaModel *GetModel() {
        return static_cast<MultiGeneSingleOmegaModel *>(model);
    }

    string GetModelType() { return modeltype; }

    //! \brief constructor for a new MCMC
    //!
    //! \param indatafile: name of file contanining sequence alignment
    //! \param intreefile: name of file contaning tree
    //! \param inevery: thinning factor
    //! \param inuntil: maximum MCMC sample size
    //! \param name: base name for all files related to this MCMC run
    //! \param force: overwrite existing files with same name
    //! \param inmyid, int innprocs: process id and total number of MPI processes
    MultiGeneSingleOmegaChain(string indatafile, string intreefile, int inblmode, int innucmode,
        int inomegamode, double inomegahypermean, double inomegahyperinvshape, int inevery,
        int inuntil, string inname, int force, int inmyid, int innprocs)
        : modeltype("MULTIGENESINGLEOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          myid(inmyid),
          nprocs(innprocs) {
        blmode = inblmode; //
        nucmode = innucmode; //
        omegamode = inomegamode; //
        omegahypermean = inomegahypermean; //
        omegahyperinvshape = inomegahyperinvshape; //
        every = inevery;                           //
        until = inuntil;                           //
        name = inname;                             //
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneSingleOmegaChain(string filename, int inmyid, int innprocs)
        : myid(inmyid), nprocs(innprocs) {
        name = filename;
        Open();
        Save();
    }

    void Reset(int force) {
        size = 0;
        if (!myid) {
            MakeFiles(force);
        }
        Save();
    }

    void New(int force) {
        model = new MultiGeneSingleOmegaModel(datafile, treefile, myid, nprocs); //
        GetModel()->SetAcrossGenesModes(blmode, nucmode, omegamode); //
        GetModel()->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape); //
        if (!myid) { cerr << "allocate\n"; } //
        GetModel()->Allocate();              //
        if (!myid) { cerr << "update\n"; }   //
        GetModel()->Update();                //
        Reset(force);
        if (!myid) { model->Trace(cerr); }
    }

    void Open() {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        is >> blmode >> nucmode >> omegamode;
        is >> omegahypermean >> omegahyperinvshape;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENESINGLEOMEGA") {
            model = new MultiGeneSingleOmegaModel(datafile, treefile, myid, nprocs);
            GetModel()->SetAcrossGenesModes(blmode, nucmode, omegamode);
            GetModel()->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape);
        } else {
            cerr << "error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        if (!myid) { cerr << "allocate\n"; }
        GetModel()->Allocate();
        model->FromStream(is);
        if (!myid) { cerr << "update\n"; }
        model->Update();
        if (!myid) {
            cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
            model->Trace(cerr);
        }
    }

    void Save() {
        if (!myid) {
            ofstream param_os((name + ".param").c_str());
            param_os << GetModelType() << '\n';
            param_os << datafile << '\t' << treefile << '\n';
            param_os << blmode << '\t' << nucmode << '\t' << omegamode << '\n';
            param_os << omegahypermean << '\t' << omegahyperinvshape << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) {
        if (myid) {
            cerr << "error: in specialized makefiles\n";
            exit(1);
        }
        if (ifstream((name + ".param").c_str()) && (force == 0)) {
            cerr << "already existing chain, cannot override (unless in forcing mode)\n";
            exit(1);
        }
        ofstream param_os((name + ".param").c_str());
        if (saveall) { ofstream chain_os((name + ".chain").c_str()); }
        ofstream mon_os((name + ".monitor").c_str());
        ofstream trace_os((name + ".trace").c_str());
        model->TraceHeader(trace_os);
        ofstream nameos((name + ".genelist").c_str());
        GetModel()->PrintGeneList(nameos);
        nameos.close();
        ofstream os((name + ".geneom").c_str());
    }

    void SavePoint() {
        if (saveall) {
            if (!myid) {
                ofstream chain_os((name + ".chain").c_str(), ios_base::app);
                GetModel()->MasterToStream(chain_os);
            } else {
                GetModel()->SlaveToStream();
            }
        }
        size++;
        if (!myid) {
            ofstream os((name + ".geneom").c_str(), ios_base::app);
            GetModel()->TraceOmega(os);
        }
    }

    void Start() {
        if (!myid) {
            ofstream run_os((name + ".run").c_str());
            run_os << 1 << '\n';
            run_os.close();
        }
        Run();
    }

    void Run() {
        if (!myid) {
            while ((GetRunningStatus() != 0) && ((until == -1) || (size <= until))) {
                MasterSendRunningStatus(1);
                Chrono chrono;
                chrono.Start();
                Move();
                chrono.Stop();

                ofstream check_os((name + ".time").c_str());
                check_os << chrono.GetTime() << '\n';
            }
            MasterSendRunningStatus(0);
            ofstream run_os((name + ".run").c_str());
            run_os << 0 << '\n';
        } else {
            while (SlaveReceiveRunningStatus()) {
                Move();
            }
        }
    }

    int GetRunningStatus() {
        ifstream run_is((name + ".run").c_str());
        int run;
        run_is >> run;
        return run;
    }
    void MasterSendRunningStatus(int status) {
        MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int SlaveReceiveRunningStatus() {
        int status;
        MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
        return status;
    }

    //! return current size (number of points saved to file thus far)
    int GetSize() { return size; }
};

int main(int argc, char *argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // int blockcounts[2] = {1, 3};
    // MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
    // MPI_Aint dtex, displacements[2];

    // displacements[0] = (MPI_Aint)0;
    // MPI_Type_extent(MPI_DOUBLE, &dtex);
    // displacements[1] = dtex;
    // MPI_Type_struct(2, blockcounts, displacements, types, &Propagate_arg);
    // MPI_Type_commit(&Propagate_arg);

    MultiGeneSingleOmegaChain *chain = 0;
    string name = "";

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneSingleOmegaChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = ""; //
        string treefile = ""; //
        int force = 1;        //
        int every = 1;        //
        int until = -1;       //
        int blmode = 1;       //
        int nucmode = 1;      //
        int omegamode = 1;    //
        double omegahypermean = 1.0; //
        double omegahyperinvshape = 1.0; //

        try {
            if (argc == 1) { throw(0); } //

            int i = 1; //
            while (i < argc) { //
                string s = argv[i]; //

                if (s == "-d") { //
                    i++;         //
                    datafile = argv[i]; //
                } else if ((s == "-t") || (s == "-T")) { //
                    i++;                                 //
                    treefile = argv[i];                  //
                } else if (s == "-f") {                  //
                    force = 1;                           //
                } else if (s == "-omega") {              //
                    omegamode = 0;                       //
                    i++;                                 //
                    string tmp = argv[i];                //
                    if (tmp != "uninf") {                //
                        omegahypermean = atof(argv[i]);  //
                        i++;                             //
                        omegahyperinvshape = atof(argv[i]); //
                    }                                       //
                } else if (s == "-nucrates") {              //
                    i++;                                    //
                    string tmp = argv[i];                   //
                    if (tmp == "shared") {                  //
                        nucmode = 2;                        //
                    } else if (tmp == "shrunken") {         //
                        nucmode = 1;                        //
                    } else if ((tmp == "ind") || (tmp == "independent")) { //
                        nucmode = 0; //
                    } else {         //
                        cerr << "error: does not recongnize command after -nucrates\n"; //
                        exit(1); //
                    }            //
                } else if (s == "-bl") { //
                    i++;                 //
                    string tmp = argv[i]; //
                    if (tmp == "shared") { //
                        blmode = 2;        //
                    } else if (tmp == "shrunken") { //
                        blmode = 1;                 //
                    } else if ((tmp == "ind") || (tmp == "independent")) { //
                        blmode = 0; //
                    } else {        //
                        cerr << "error: does not recongnize command after -bl\n"; //
                        exit(1); //
                    }            //
                } else if ((s == "-x") || (s == "-extract")) { //
                    i++;                                       //
                    if (i == argc) throw(0);                   //
                    every = atoi(argv[i]);                     //
                    i++;                                       //
                    if (i == argc) throw(0);                   //
                    until = atoi(argv[i]);                     //
                } else {                                       //
                    if (i != (argc - 1)) { throw(0); }         //
                    name = argv[i];                            //
                }                                              //
                i++;                                           //
            }                                                  //
            if ((datafile == "") || (treefile == "") || (name == "")) { throw(0); } //
        } catch (...) { //
            cerr << "globom -d <alignment> -t <tree> <chainname> \n"; //
            cerr << '\n'; //
            exit(1);      //
        }                 //

        chain = new MultiGeneSingleOmegaChain(datafile, treefile, blmode, nucmode, omegamode,
            omegahypermean, omegahyperinvshape, every, until, name, force, myid, nprocs);
    }

    if (!myid) { cerr << "chain " << name << " started\n"; }
    chain->Start();
    if (!myid) {
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }

    MPI_Finalize();
}