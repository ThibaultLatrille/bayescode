
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "IIDDirichlet.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"
#include "StickBreakingProcess.hpp"
#include "MultinomialAllocationVector.hpp"

const int Nrr = Nnuc * (Nnuc-1) / 2;
const int Nstate = 61;

class AAMutSelHyperSBDPOmegaModel : public ProbModel {

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

	std::vector<double> nucstat;
	std::vector<double> nucrelrate;
	GTRSubMatrix* nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
	double omega;
	OmegaSuffStat omegasuffstat;
	
    // mixture components
    // set of Ncat Dirichlet densities
    // centers
    vector<double> aacenterhypercenter;
    double aacenterhyperinvconc;
    IIDDirichlet* componentaacenterarray;
    // concentrations
    double aaconchypermean;
    double aaconchyperinvshape;
    IIDGamma* componentaaconcentrationarray;
    // and associated suffstatarray
    DirichletSuffStatArray* aahypersuffstatarray;

    // now, a mixture model drawing from this set of Ncat components
    // weights:
    double kappa;
    StickBreakingProcess* weight;

    // site allocations
	MultinomialAllocationVector* sitealloc;
    // dispatching across sites
    ConstMixtureArray<vector<double> >* siteaacenterarray;
    ConstMixtureArray<double>* siteaaconcentrationarray;

    // aa fitness arrays:w
    //
    MultiDirichlet* aafitnessarray;

	AAMutSelOmegaCodonSubMatrixArray* codonmatrixarray;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* pathsuffstatarray;

    int Ncat;

    int fixhypermix;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	AAMutSelHyperSBDPOmegaModel(string datafile, string treefile, int inNcat) {

        fixhypermix = 0;

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

        Ncat = inNcat;
        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100)    {
                Ncat = 100;
            }
        }


		std::cerr << "-- Number of sites: " << Nsite << std::endl;

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		tree->SetIndices();
		Nbranch = tree->GetNbranch();

		// Allocate();
	}

    void SetFixAAHyperMix(int inmix)    {
        fixhypermix = inmix;
    }

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogProb() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
		// Trace(cerr);
    }

	void Allocate()	{

		lambda = 10;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);
		lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);

		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));

		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        // mixture components
        aacenterhypercenter.assign(Naa,1.0/Naa);
        aacenterhyperinvconc = 1.0/Naa;
        componentaacenterarray = new IIDDirichlet(Ncat,aacenterhypercenter,1.0/aacenterhyperinvconc);
        componentaacenterarray->SetUniform();

        aaconchypermean = Naa;
        aaconchyperinvshape = 1.0;
        double alpha = 1.0 / aaconchyperinvshape;
        double beta = alpha / aaconchypermean;
        componentaaconcentrationarray = new IIDGamma(Ncat,alpha,beta);
        for (int k=0; k<Ncat; k++)  {
            (*componentaaconcentrationarray)[k] = 20.0;
        }

        // suff stats for mixture components
        aahypersuffstatarray = new DirichletSuffStatArray(Ncat,Naa);

        // mixture weights
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat,kappa);

        // site allocations
        sitealloc = new MultinomialAllocationVector(Nsite,weight->GetArray());

        // site distributors
        siteaacenterarray = new ConstMixtureArray<vector<double> >(componentaacenterarray,sitealloc);
        siteaaconcentrationarray = new ConstMixtureArray<double>(componentaaconcentrationarray,sitealloc);

        // site-specific amino-acid fitness profiles
        aafitnessarray = new MultiDirichlet(siteaacenterarray,siteaaconcentrationarray);

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        codonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, aafitnessarray, omega);
		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,codonmatrixarray);
		pathsuffstatarray = new PathSuffStatArray(Nsite);
	}

    //-------------------
    // Accessors
    // ------------------

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    double GetOmega() const {
        return omega;
    }

    const PoissonSuffStatBranchArray* GetLengthSuffStatArray() const {
        return lengthsuffstatarray;
    }

    const DirichletSuffStatArray* GetAAHyperSuffStatArray() const   {
        return aahypersuffstatarray;
    }

    const vector<int>& GetOccupancies() const   {
        return sitealloc->GetOccupancies();
    }

    //-------------------
    // Setting and updating
    // ------------------

    void SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    void SetOmega(double inomega)   {
        omega = inomega;
        UpdateCodonMatrices();
    }

    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape)   {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    void SetAAHyperMixture(const ConstArray<vector<double> >& incomponentaacenterarray, const ConstArray<double>& incomponentaaconcentrationarray, const ConstArray<double>& inweight) {

        componentaacenterarray->Copy(incomponentaacenterarray);
        componentaaconcentrationarray->Copy(incomponentaaconcentrationarray);
        weight->Copy(inweight);
    }

    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

	void UpdateCodonMatrices()	{
        codonmatrixarray->SetOmega(omega);
		codonmatrixarray->UpdateCodonMatrices();
	}

    void UpdateCodonMatrix(int site)    {
        (*codonmatrixarray)[site].CorruptMatrix();
    }
		
    void UpdateMatrices()   {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    void NoUpdate() {}

    //-------------------
    // Priors and likelihood
    //-------------------

    double GetLogPrior() const {
        double total = 0;
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AAHyperLogPrior();
        total += AALogPrior();
        total += OmegaLogPrior();
        return total;
    }

	double GetLogLikelihood() const {
		return phyloprocess->GetLogProb();
	}

    double GetLogProb() const   {
        return GetLogPrior() + GetLogLikelihood();
    }

	double BranchLengthsHyperLogPrior() const {
		return -lambda / 10;
	}

	double BranchLengthsLogPrior() const {
		return branchlength->GetLogProb();
	}

    double StickBreakingHyperLogPrior() const   {
        return -kappa/10;
    }

    double StickBreakingLogPrior() const    {
        return weight->GetLogProb(kappa);
        // return weight->GetMarginalLogProb(sitealloc->GetOccupancies());
    }

	// exponential of mean 1
	double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		return alpha * log(beta) - Random::logGamma(alpha) + (alpha-1) * log(omega) - beta*omega;
	}

    double AAHyperLogPrior() const {
        double total = 0;
        total += componentaacenterarray->GetLogProb();
        total += componentaaconcentrationarray->GetLogProb();
        return total;
    }

    double AAHyperLogPrior(int k) const {
        double total = 0;
        total += componentaacenterarray->GetLogProb(k);
        total += componentaaconcentrationarray->GetLogProb(k);
        return total;
    }

    double AALogPrior() const {
        return aafitnessarray->GetLogProb();
    }

    double AALogPrior(int i) const {
        return aafitnessarray->GetLogProb(i);
    }

    double NucRatesLogPrior() const {
        return 0;
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

	double PathSuffStatLogProb() const {
		return pathsuffstatarray->GetLogProb(*codonmatrixarray);
	}

    double PathSuffStatLogProb(int site) const {
        return pathsuffstatarray->GetVal(site).GetLogProb(codonmatrixarray->GetVal(site));
    }

	double BranchLengthsHyperSuffStatLogProb() const {
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double AAHyperSuffStatLogProb(int k) const   {
        return aahypersuffstatarray->GetVal(k).GetLogProb(componentaacenterarray->GetVal(k),componentaaconcentrationarray->GetVal(k));
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // for moving branch lengths hyperparameter lambda
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // for moving nuc rates
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + PathSuffStatLogProb();
    }

    // for moving aa hyper params (aacenter and aainvconc)
    // for component k of the mixture
    double AAHyperLogProb(int k) const   {
        return AAHyperLogPrior(k) + AAHyperSuffStatLogProb(k);
    }

    // for moving kappa
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Moves 
    //-------------------

	double Move()	{
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    void ResampleSub(double frac)   {
        UpdateMatrices();
		phyloprocess->Move(frac);
    }

    void MoveParameters(int nrep)   {
		for (int rep=0; rep<nrep; rep++)	{

			ResampleBranchLengths();
			MoveBranchLengthsHyperParameter();

			CollectPathSuffStat();

            MoveAA(3);

			MoveNucRates();
			MoveOmega();
		}
	}

    void CollectLengthSuffStat()    {
		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
    }

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
        ScalingMove(lambda,1.0,10,&AAMutSelHyperSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelHyperSBDPOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelHyperSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelHyperSBDPOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		pathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(*pathsuffstatarray);
	}

	void MoveOmega()	{

		omegasuffstat.Clear();
		omegasuffstat.AddSuffStat(*codonmatrixarray,*pathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegasuffstat.GetCount(), beta + omegasuffstat.GetBeta());
		UpdateCodonMatrices();
	}

	void MoveNucRates()	{

        UpdateMatrices();

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelHyperSBDPOmegaModel::NucRatesLogProb,&AAMutSelHyperSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelHyperSBDPOmegaModel::NucRatesLogProb,&AAMutSelHyperSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelHyperSBDPOmegaModel::NucRatesLogProb,&AAMutSelHyperSBDPOmegaModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelHyperSBDPOmegaModel::NucRatesLogProb,&AAMutSelHyperSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelHyperSBDPOmegaModel::NucRatesLogProb,&AAMutSelHyperSBDPOmegaModel::UpdateMatrices,this);

        UpdateMatrices();
	}

    void MoveAA(int nrep)   {

        for (int rep=0; rep<nrep; rep++)    {
            MoveAAFitness();
            ResampleAlloc();
            if (! fixhypermix)  {
                MoveAAHyperMixture();
            }
        }
    }

    void MoveAAHyperMixture()   {
        for (int rep=0; rep<10; rep++)  {
            MoveAAHyperMixtureComponents();
            LabelSwitchingMove();
            ResampleWeights();
            MoveKappa();
        }
    }

    void MoveAAHyperMixtureComponents() {

        CollectAAHyperSuffStat();

        for (int k=0; k<Ncat; k++)  {
            for (int rep=0; rep<3; rep++)    {
                MoveAAHyperCenters(1.0,1);
                MoveAAHyperCenters(0.3,3);
                MoveAAHyperConcentrations(1.0);
                MoveAAHyperConcentrations(0.3);
            }
        }
    }

    void CollectAAHyperSuffStat()   {
        aahypersuffstatarray->Clear();
        aafitnessarray->AddSuffStat(*aahypersuffstatarray,*sitealloc);
    }

    double MoveAAHyperCenters(double tuning, int n) {
		double nacc = 0;
		double ntot = 0;
        vector<double> bk(Naa,0);
        for (int k=0; k<Ncat; k++)  {
            vector<double>& aa = (*componentaacenterarray)[k];
            bk = aa;
            double deltalogprob = -AAHyperLogProb(k);
            double loghastings = Random::ProfileProposeMove(aa,Naa,tuning,n);
            deltalogprob += loghastings;
            deltalogprob += AAHyperLogProb(k);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                aa = bk;
            }
            ntot++;
        }
		return nacc/ntot;
	}

    double MoveAAHyperConcentrations(double tuning)  {
		double nacc = 0;
		double ntot = 0;
        for (int k=0; k<Ncat; k++)  {
            double& c = (*componentaaconcentrationarray)[k];
            double bk = c;
            double deltalogprob = -AAHyperLogProb(k);
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            c *= e;
            deltalogprob += m;
            deltalogprob += AAHyperLogProb(k);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                c = bk;
            }
            ntot++;
        }
		return nacc/ntot;
    }

    void ResampleAlloc()    {
        vector<double> postprob(Ncat,0);
        for (int i=0; i<Nsite; i++) {
            GetAllocPostProb(i,postprob);
            sitealloc->GibbsResample(i,postprob);
        }
        sitealloc->UpdateOccupancies();
    }

    void GetAllocPostProb(int site, vector<double>& postprob)    {

        double max = 0;
        const vector<double>& w = weight->GetArray();
        for (int i=0; i<Ncat; i++) {
            double tmp = Random::logDirichletDensity(aafitnessarray->GetVal(site),componentaacenterarray->GetVal(i),componentaaconcentrationarray->GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp))    {
                max = tmp;
            }
        }

        double total = 0;
        for (int i=0; i<Ncat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i=0; i<Ncat; i++) {
            postprob[i] /= total;
        }
    }

    void LabelSwitchingMove()   {
        MoveOccupiedCompAlloc(5);
        MoveAdjacentCompAlloc(5);
    }

    double MoveOccupiedCompAlloc(int k0)	{

        const vector<int>& occupancy = sitealloc->GetOccupancies();
        const vector<double>& V = weight->GetBetaVariates();
        const vector<double>& w = weight->GetArray();

        int nrep = (int) (k0 * kappa);
        ResampleWeights();
        double total = 0.0;
        int Nocc = GetNcluster();
        if (Nocc != 1)	{
            for (int i=0; i<nrep; i++)	{
                int occupiedComponentIndices[Nocc];
                int j=0;
                for (int k=0; k<Ncat; k++)	{
                    if (occupancy[k] != 0)	{
                        occupiedComponentIndices[j] = k;
                        j++;
                    }
                }
                if (j != Nocc)	{
                    cerr << "error in MoveOccupiedCompAlloc.\n";
                    exit(1);
                }
                int indices[2];
                Random::DrawFromUrn(indices,2,Nocc);
                int cat1 = occupiedComponentIndices[indices[0]];
                int cat2 = occupiedComponentIndices[indices[1]];
                double logMetropolis = (occupancy[cat2] - occupancy[cat1]) * log(w[cat1] / w[cat2]);
                int accepted = (log(Random::Uniform()) < logMetropolis);
                if (accepted)	{
                    total += 1.0;
                    componentaacenterarray->Swap(cat1,cat2);
                    componentaaconcentrationarray->Swap(cat1,cat2);
                    sitealloc->SwapComponents(cat1,cat2);
                }
            }
            return total /= nrep;
        }
        return 0;
    }

    double MoveAdjacentCompAlloc(int k0)	{

        ResampleWeights();
        int nrep = (int) (k0 * kappa);
        
        double total = 0;

        const vector<int>& occupancy = sitealloc->GetOccupancies();
        const vector<double>& V = weight->GetBetaVariates();

        for (int i=0; i<nrep; i++)	{
            int cat1 = (int)(Random::Uniform() * (Ncat-2));  
            int cat2 = cat1 + 1;
            double logMetropolis = (occupancy[cat1] * log(1 - V[cat2])) - (occupancy[cat2] * log(1-V[cat1]));
            int accepted = (log(Random::Uniform()) < logMetropolis);
            if (accepted)	{
                total += 1.0;
                componentaacenterarray->Swap(cat1,cat2);
                componentaaconcentrationarray->Swap(cat1,cat2);
                sitealloc->SwapComponents(cat1,cat2);
                weight->SwapComponents(cat1,cat2);
            }
        }

        return total /= nrep;
    }

    void ResampleWeights()  {
        weight->GibbsResample(sitealloc->GetOccupancies());
    }

    void MoveKappa()    {
        ScalingMove(kappa,1.0,10,&AAMutSelHyperSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelHyperSBDPOmegaModel::NoUpdate,this);
        ScalingMove(kappa,0.3,10,&AAMutSelHyperSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelHyperSBDPOmegaModel::NoUpdate,this);
    }

    int GetNcluster() const {

        int n = 0;
        const vector<int>& occupancy = sitealloc->GetOccupancies();
        for (int i=0; i<Ncat; i++)  {
            if (occupancy[i])    {
                n++;
            }
        }
        return n;
    }

    double MoveAAFitness() {
        MoveAAFitness(1.0,1,3);
        MoveAAFitness(0.3,1,3);
        MoveAAFitness(0.1,3,3);
        MoveAAFitness(0.1,5,3);
        return 1.0;
    }

	double MoveAAFitness(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Naa];
        for (int i=0; i<Nsite; i++) {
            vector<double>& aa = (*aafitnessarray)[i];
            for (int rep=0; rep<nrep; rep++)	{
                for (int l=0; l<Naa; l++)	{
                    bk[l] = aa[l];
                }
                double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                double loghastings = Random::ProfileProposeMove(aa,Naa,tuning,n);
                deltalogprob += loghastings;
                UpdateCodonMatrix(i);
                deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted)	{
                    nacc ++;
                }
                else	{
                    for (int l=0; l<Naa; l++)	{
                        aa[l] = bk[l];
                    }
                    UpdateCodonMatrix(i);
                }
                ntot++;
            }
        }
		return nacc/ntot;
	}

    //-------------------
    // Traces and Monitors
    // ------------------

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "aaent\t";
        os << "ncluster\t";
        os << "kappa\t";
        os << "aahypercenter\t";
		os << "hyperstatent\t";
		os << "rrent\t";
        os << "statent\n";
	}

	void Trace(ostream& os) const {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
		os << omega << '\t';
        os << aafitnessarray->GetMeanEntropy() << '\t';
        os << GetNcluster() << '\t';
        os << kappa << '\t';
        os << componentaacenterarray->GetMeanEntropy() << '\t';
		os << componentaaconcentrationarray->GetMean() << '\t';
		os << Random::GetEntropy(nucrelrate) << '\t';
		os << Random::GetEntropy(nucstat) << '\n';
	}

	void Monitor(ostream& os) const {}

	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

};

