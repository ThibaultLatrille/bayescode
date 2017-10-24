
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

class AAMutSelSBDPOmegaModel : public ProbModel {

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

    int Ncat;

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
	
    // mixture over amino-acid fitness profiles

    // (1) component weights (truncated stick breaking)
    double kappa;
    StickBreakingProcess* weight;

    // (2) component values: independent aa fitness profiles from Dirichlet distribution
    vector<double> aacenter;
    double aainvconc;
    IIDDirichlet* componentaafitnessarray;
    DirichletSuffStat aahypersuffstat;

    // (3) site allocations: multinomial given the weights
    // multinomial allocation of sites to components of aa fitness profile distribution
	MultinomialAllocationVector* sitealloc;

    // an array of codon matrices (one for each distinct aa fitness profile)
	AAMutSelOmegaCodonSubMatrixArray* componentcodonmatrixarray;

	// this one is used by PhyloProcess: has to be a ConstArray<SubMatrix>
	ConstMixtureArray<SubMatrix>* sitesubmatrixarray;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* sitepathsuffstatarray;
	PathSuffStatArray* componentpathsuffstatarray;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	AAMutSelSBDPOmegaModel(string datafile, string treefile, int inNcat) : aahypersuffstat(20) {

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

        Ncat = inNcat;

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

        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat,kappa);

        sitealloc = new MultinomialAllocationVector(Nsite,weight->GetArray());

        aacenter.assign(Naa,1.0/Naa);
        aainvconc = 1.0/Naa;
        componentaafitnessarray = new IIDDirichlet(Ncat,aacenter,1.0/aainvconc);

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        componentcodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, componentaafitnessarray, omega);

        sitesubmatrixarray = new ConstMixtureArray<SubMatrix>(componentcodonmatrixarray,sitealloc);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitesubmatrixarray);
		sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatarray = new PathSuffStatArray(Ncat);
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
        componentcodonmatrixarray->SetOmega(omega);
		componentcodonmatrixarray->UpdateCodonMatrices();
	}

    void UpdateCodonMatrix(int k)    {
        (*componentcodonmatrixarray)[k].CorruptMatrix();
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

	// exponential of mean 1
	double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		return alpha * log(beta) - Random::logGamma(alpha) + (alpha-1) * log(omega) - beta*omega;
	}

    double StickBreakingHyperLogPrior() const   {
        return -kappa/10;
    }

    double StickBreakingLogPrior() const    {
        return weight->GetLogProb(kappa);
        // return weight->GetMarginalLogProb(sitealloc->GetOccupancies());
    }

    double AAHyperLogPrior() const {
        return -aainvconc;
    }

    double AALogPrior() const {
        return componentaafitnessarray->GetLogProb();
    }

    double AALogPrior(int k) const {
        return componentaafitnessarray->GetLogProb(k);
    }

    double NucRatesLogPrior() const {
        return 0;
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

	double PathSuffStatLogProb() const {
        return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
	}

    double PathSuffStatLogProb(int k) const {
        return componentpathsuffstatarray->GetVal(k).GetLogProb(componentcodonmatrixarray->GetVal(k));
    }

	double BranchLengthsHyperSuffStatLogProb() const {
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double AAHyperSuffStatLogProb() const   {
        return aahypersuffstat.GetLogProb(aacenter,1.0/aainvconc);
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
    double AAHyperLogProb() const   {
        return AAHyperLogPrior() + AAHyperSuffStatLogProb();
    }

    // for moving kappa
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    // per site
	void CollectSitePathSuffStat()	{
		sitepathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
	}

    // per component of the mixture
	void CollectComponentPathSuffStat()	{
        componentpathsuffstatarray->Clear();
        sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
    }

    void CollectLengthSuffStat()    {
		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
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

			CollectSitePathSuffStat();
            CollectComponentPathSuffStat();

			MoveNucRates();

			MoveOmega();

            MoveAAMixture();
            MoveAAHyperParameters();
            MoveKappa();
		}
	}

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{
		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
        ScalingMove(lambda,1.0,10,&AAMutSelSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelSBDPOmegaModel::BranchLengthsHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void MoveOmega()	{

		omegasuffstat.Clear();
		omegasuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegasuffstat.GetCount(), beta + omegasuffstat.GetBeta());
		UpdateCodonMatrices();
	}

	void MoveNucRates()	{

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelSBDPOmegaModel::NucRatesLogProb,&AAMutSelSBDPOmegaModel::UpdateMatrices,this);
	}

    void MoveAAHyperParameters()    {

        aahypersuffstat.Clear();
        componentaafitnessarray->AddSuffStat(aahypersuffstat,sitealloc->GetOccupancies());
        for (int rep=0; rep<10; rep++)  {
            ProfileMove(aacenter,0.1,1,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ProfileMove(aacenter,0.03,3,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ProfileMove(aacenter,0.01,3,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ScalingMove(aainvconc,1.0,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
            ScalingMove(aainvconc,0.3,10,&AAMutSelSBDPOmegaModel::AAHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
        }
        componentaafitnessarray->SetCenter(aacenter);
        componentaafitnessarray->SetConcentration(1.0/aainvconc);
        componentaafitnessarray->PriorResample(sitealloc->GetOccupancies());
        componentcodonmatrixarray->UpdateCodonMatrices(sitealloc->GetOccupancies());
    }

    void MoveAAMixture()    {

        for (int rep=0; rep<10; rep++)  {
            MoveAAProfiles();
            ResampleAlloc();
            LabelSwitchingMove();
            ResampleWeights();
            UpdateCodonMatrices();
            CollectComponentPathSuffStat();
        }
    }

    double MoveAAProfiles() {
        MoveAAProfiles(1.0,1,3);
        MoveAAProfiles(0.3,1,3);
        MoveAAProfiles(0.1,3,3);
        MoveAAProfiles(0.1,5,3);
        componentaafitnessarray->PriorResample(sitealloc->GetOccupancies());
        componentcodonmatrixarray->UpdateCodonMatrices(sitealloc->GetOccupancies());
        return 1.0;
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
        const PathSuffStat& suffstat = sitepathsuffstatarray->GetVal(site);
        for (int i=0; i<Ncat; i++) {
            double tmp = suffstat.GetLogProb(componentcodonmatrixarray->GetVal(i));
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

    void SwapComponents(int cat1, int cat2) {
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
                    componentaafitnessarray->Swap(cat1,cat2);
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
                componentaafitnessarray->Swap(cat1,cat2);
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
        ScalingMove(kappa,1.0,10,&AAMutSelSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(kappa,0.3,10,&AAMutSelSBDPOmegaModel::StickBreakingHyperLogProb,&AAMutSelSBDPOmegaModel::NoUpdate,this);
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

	double MoveAAProfiles(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Naa];
        const vector<int>& occupancy = sitealloc->GetOccupancies();
        for (int i=0; i<Ncat; i++) {
            if (occupancy[i])   {
                vector<double>& aa = (*componentaafitnessarray)[i];
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
        }
		return nacc/ntot;
	}

    //-------------------
    // Traces and Monitors
    // ------------------

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "ncluster\t";
        os << "kappa\t";
        os << "aaent\t";
        os << "aainvconc\t";
        os << "aacenter\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) const {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
		os << omega << '\t';
        os << GetNcluster() << '\t';
        os << kappa << '\t';
        os << componentaafitnessarray->GetMeanEntropy() << '\t';
        os << aainvconc << '\t';
        os << Random::GetEntropy(aacenter) << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
	}

	void Monitor(ostream& os) const {}

	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

};

