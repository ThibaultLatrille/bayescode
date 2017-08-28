#include "CodonM2aModel.hpp"

// constuction allocation
//

CodonM2aModel::CodonM2aModel(string datafile, string treefile, double inpi)	{

    blmode = 0;
    nucmode = 0;
    data = new FileSequenceAlignment(datafile);
    codondata = new CodonSequenceAlignment(data, true);
    pi = inpi;

    Nsite = codondata->GetNsite();    // # columns
    Ntaxa = codondata->GetNtaxa();

    // std::cerr << "-- Number of sites: " << Nsite << std::endl;

    taxonset = codondata->GetTaxonSet();

    // get tree from file (newick format)
    tree = new Tree(treefile);

    // check whether tree and data fits together
    tree->RegisterWith(taxonset);

    tree->SetIndices();
    Nbranch = tree->GetNbranch();
}


void CodonM2aModel::Unfold()   {

    phyloprocess->Unfold();
    phyloprocess->ResampleSub();
}

void CodonM2aModel::Allocate()	{

    lambda = 10.0;
    blhypermean = new BranchIIDGamma(*tree,1.0,lambda);
    blhypermean->SetAllBranches(1.0 / lambda);
    blhyperinvshape = 1.0;
    branchlength = new GammaWhiteNoise(*tree,*blhypermean,1.0/blhyperinvshape);

    purom = 0.5;
    puromhypermean = 0.5;
    puromhyperinvconc = 0.5;

    dposom = 1.0;
    dposomhypermean = 0.5;
    dposomhyperinvshape = 0.5;

    purw = 0.1;
    purwhypermean = 0.5;
    purwhyperinvconc = 0.5;

    if (! pi)   {
        posw = 0;
        poswhypermean = 0;
        poswhyperinvconc = 0;
    }
    else    {
        posw = 0.1;
        poswhypermean = 0.5;
        poswhyperinvconc = 0.5;
    }

    componentomegaarray = new M2aMix(purom,dposom+1,purw,posw);
    sitealloc = new MultinomialAllocationVector(GetNsite(),componentomegaarray->GetWeights());
    sitepostprobarray.assign(GetNsite(),vector<double>(3,0));

    nucrelratehypercenter.assign(Nrr,1.0/Nrr);
    nucrelratehyperinvconc = 1.0 / Nrr;

    nucstathypercenter.assign(Nnuc,1.0/Nnuc);
    nucstathyperinvconc = 1.0 / Nnuc;

    nucrelrate.assign(Nrr,0);
    double totrr = 0;
    for (int k=0; k<Nrr; k++)	{
        nucrelrate[k] = Random::sExpo();
        totrr += nucrelrate[k];
    }
    for (int k=0; k<Nrr; k++)	{
        nucrelrate[k] /= totrr;
    }

    nucstat.assign(Nnuc,0);
    double totstat = 0;
    for (int k=0; k<Nnuc; k++)	{
        nucstat[k] = Random::sGamma(1.0);
        totstat += nucstat[k];
    }
    for (int k=0; k<Nnuc; k++)	{
        nucstat[k] /= totstat;
    }
    nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

    componentcodonmatrixarray = new MGOmegaCodonSubMatrixArray((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,componentomegaarray);

    sitesubmatrixarray = new ConstMixtureArray<SubMatrix>(componentcodonmatrixarray,sitealloc);
    sitecodonmatrixarray = new ConstMixtureArray<MGOmegaCodonSubMatrix>(componentcodonmatrixarray,sitealloc);

    phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitesubmatrixarray);

    lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);
    sitepathsuffstatarray = new PathSuffStatArray(GetNsite());
    componentpathsuffstatarray = new PathSuffStatArray(3);
    siteomegasuffstatarray = new OmegaSuffStatArray(GetNsite());
}

// setting model features and (hyper)parameters
//

void CodonM2aModel::SetAcrossGenesModes(int inblmode, int innucmode)   {
    blmode = inblmode;
    nucmode = innucmode;
}

void CodonM2aModel::SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
    for (int j=0; j<Nbranch; j++)   {
        (*branchlength)[j] = inbranchlength.GetVal(j);
    }
}

void CodonM2aModel::GetBranchLengths(BranchArray<double>& inbranchlength) const   {
    for (int j=0; j<Nbranch; j++)   {
        inbranchlength[j] = branchlength->GetVal(j);
    }
}

void CodonM2aModel::SetBranchLengthsHyperParameters(const ConstBranchArray<double>& inblmean, double inblinvshape) {
    for (int j=0; j<Nbranch; j++)   {
        (*blhypermean)[j] = inblmean.GetVal(j);
    }
    blhyperinvshape = inblinvshape;
    branchlength->SetShape(1.0 / blhyperinvshape);
}

void CodonM2aModel::SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {

    nucrelrate = innucrelrate;
    nucstat = innucstat;
    UpdateMatrices();
}

void CodonM2aModel::GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const {

    innucrelrate = nucrelrate;
    innucstat = nucstat;
}

void CodonM2aModel::SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc, const std::vector<double>& innucstathypercenter, double innucstathyperinvconc) {

    nucrelratehypercenter = innucrelratehypercenter;
    nucrelratehyperinvconc = innucrelratehyperinvconc;
    nucstathypercenter = innucstathypercenter;
    nucstathyperinvconc = innucstathyperinvconc;
}

void CodonM2aModel::SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw)    {

    purom = inpurom;
    dposom = indposom;
    purw = inpurw;
    posw = inposw;
    componentomegaarray->SetParameters(purom,dposom+1,purw,posw);
}

void CodonM2aModel::GetMixtureParameters(double& inpurom, double& indposom, double& inpurw, double& inposw)    {

    inpurom = purom;
    indposom = dposom;
    inpurw = purw;
    inposw = posw;
}

void CodonM2aModel::SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape, double inpi, double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc)  {

    puromhypermean = inpuromhypermean;
    puromhyperinvconc = inpuromhyperinvconc;
    dposomhypermean = indposomhypermean;
    dposomhyperinvshape = indposomhyperinvshape;
    pi = inpi;
    purwhypermean = inpurwhypermean;
    purwhyperinvconc = inpurwhyperinvconc;
    poswhypermean = inposwhypermean;
    poswhyperinvconc = inposwhyperinvconc;
}

// 
// Matrices
//

void CodonM2aModel::UpdateNucMatrix()	{
    nucmatrix->CopyStationary(nucstat);
    nucmatrix->CorruptMatrix();
}

void CodonM2aModel::UpdateCodonMatrices()	{
    componentcodonmatrixarray->UpdateCodonMatrices();
}
    
void CodonM2aModel::UpdateMatrices()   {
    UpdateNucMatrix();
    UpdateCodonMatrices();
}

//
// Likelihood
//

double CodonM2aModel::GetLogLikelihood()	{
    // return GetIntegratedLogLikelihood();
    return phyloprocess->GetLogProb();
}

double CodonM2aModel::GetIntegratedLogLikelihood() {

    int ncat = 3;

    double total = 0;
    double logp[ncat];
    const vector<double>& w = componentomegaarray->GetWeights();
    double max = 0;
    for (int i=0; i<GetNsite(); i++) {
        int bkalloc = sitealloc->GetVal(i);

        for (int k=0; k<ncat; k++) {
            (*sitealloc)[i] = k;
            logp[k] = phyloprocess->SiteLogLikelihood(i);
            if ((!k) || (max<logp[k]))  {
                max = logp[k];
            }
        }

        double p = 0;
        for (int k=0; k<ncat; k++) {
            p += w[k] * exp(logp[k]-max);
        }
        double logl = log(p) + max;
        total += logl;

        (*sitealloc)[i] = bkalloc;
    }
    return total;
}

//
// Suff Stat and suffstatlogprobs
//

const PoissonSuffStatBranchArray* CodonM2aModel::GetLengthSuffStatArray()  {
    return lengthsuffstatarray;
}

double CodonM2aModel::LambdaHyperSuffStatLogProb()	{
    return lambdasuffstat.GetLogProb(1.0,lambda);
}

const NucPathSuffStat& CodonM2aModel::GetNucPathSuffStat() {
    return nucpathsuffstat;
}

double CodonM2aModel::NucRatesSuffStatLogProb() {
    return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
}

double CodonM2aModel::PathSuffStatLogProb()	{
    return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
}

double CodonM2aModel::OmegaSuffStatLogProb()    {
    componentomegaarray->SetParameters(purom,dposom+1,purw,posw);
    return componentomegaarray->GetPostProbArray(*siteomegasuffstatarray,sitepostprobarray);
}

//
// Priors
//

double CodonM2aModel::GetLogPrior() {
    double total = 0;
    total += BranchLengthsLogPrior();
    total += NucRatesLogPrior();
    total += OmegaLogPrior();
    return total;
}

double CodonM2aModel::LambdaHyperLogPrior()	{
    return -lambda / 10;
}

double CodonM2aModel::BranchLengthsLogPrior()   {

    double total = 0;
    if (blmode == 0)    {
        total += LambdaHyperLogPrior();
    }
    total += branchlength->GetLogProb();
    return total;
}

double CodonM2aModel::NucRatesLogPrior()    {

    double total = 0;
    double rrconc = 1.0 / nucrelratehyperinvconc;
    total += Random::logGamma(rrconc);
    for (int i=0; i<Nrr; i++)   {
        total += -Random::logGamma(rrconc * nucrelratehypercenter[i]) + (rrconc*nucrelratehypercenter[i] -1)*log(nucrelrate[i]);
    }
    double statconc = 1.0 / nucstathyperinvconc;
    total += Random::logGamma(statconc);
    for (int i=0; i<Nnuc; i++)  {
        total += -Random::logGamma(statconc * nucstathypercenter[i]) + (statconc*nucstathypercenter[i]-1)*log(nucstat[i]);
    }
    return total;
}

//
// Hyper priors for omega mixture
//

double CodonM2aModel::OmegaLogPrior()  {
    double total = 0;
    total += PurOmegaLogProb();
    total += PosOmegaLogProb();
    total += PurWeightLogProb();
    total += PosWeightLogProb();
    return total;
}

// Beta prior for purifmean
double CodonM2aModel::PurOmegaLogProb()   {
    double alpha = puromhypermean / puromhyperinvconc;
    double beta = (1-puromhypermean) / puromhyperinvconc;
    return Random::logGamma(alpha+beta) - Random::logGamma(alpha) - Random::logGamma(beta) + (alpha-1)*log(purom) + (beta-1)*log(1.0-purom);
}

// Gamma prior for dposom
double CodonM2aModel::PosOmegaLogProb()	{
    double alpha = 1.0 / dposomhyperinvshape;
    double beta = alpha / dposomhypermean;
    return alpha*log(beta) - Random::logGamma(alpha) + (alpha-1)*log(dposom) - beta*dposom;
}

// Beta prior for purw
double CodonM2aModel::PurWeightLogProb()   {
    double alpha = purwhypermean / purwhyperinvconc;
    double beta = (1 - purwhypermean) / purwhyperinvconc;
    return Random::logGamma(alpha+beta) - Random::logGamma(alpha) - Random::logGamma(beta) + (alpha-1)*log(purw) + (beta-1)*log(1.0-purw);
}

// mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1 - pi) for posw
double CodonM2aModel::PosWeightLogProb()   {
    if (posw)   {
        if (! pi)   {
            cerr << "in PosWeightLogProb: pi == 0 and posw > 0\n";
            exit(1);
        }

        double alpha = poswhypermean / poswhyperinvconc;
        double beta = (1 - poswhypermean) / poswhyperinvconc;
        return log(pi) + Random::logGamma(alpha+beta) - Random::logGamma(alpha) - Random::logGamma(beta) + (alpha-1)*log(posw) + (beta-1)*log(1.0-posw);
    }
    else    {
        return log(1-pi);
    }
}

// Bernoulli for whether posw == 0 or > 0
double CodonM2aModel::PosSwitchLogProb()   {
    if (posw)   {
        return log(pi);
    }
    return log(1-pi);
}

//
//  Moves 
//

void CodonM2aModel::Move()	{

    Chrono mappingtime, reptime, lengthtime, collecttime, omegatime, nuctime;

    mappingtime.Start();
    ResampleSub(1.0);
    mappingtime.Stop();

    int nrep = 30;

    reptime.Start();
    for (int rep=0; rep<nrep; rep++)	{

        if (blmode != 2)    {
            lengthtime.Start();
            MoveBranchLengths();
            lengthtime.Stop();
        }

        collecttime.Start();
        CollectPathSuffStat();
        collecttime.Stop();

        omegatime.Start();
        MoveOmega();
        omegatime.Stop();

        if (nucmode != 2)    {
            nuctime.Start();
            UpdateMatrices();
            MoveNucRates();
            nuctime.Stop();
        }
    }
    reptime.Stop();

    // cerr << mappingtime.GetTime() << '\t' << reptime.GetTime() << '\t' << lengthtime.GetTime() + collecttime.GetTime() + omegatime.GetTime() + nuctime.GetTime() << '\t' << lengthtime.GetTime() << '\t' << collecttime.GetTime() << '\t' << omegatime.GetTime() << '\t' << nuctime.GetTime() << '\n';
}

void CodonM2aModel::ResampleSub(double frac)  {
    UpdateMatrices();
    phyloprocess->Move(frac);
}

//
// Branch Lengths and hyperparam lambda
//

void CodonM2aModel::MoveBranchLengths()    {
        ResampleBranchLengths();
        if (blmode == 0)    {
            MoveLambda();
        }
}

void CodonM2aModel::ResampleBranchLengths()	{

    CollectLengthSuffStat();
    branchlength->GibbsResample(*lengthsuffstatarray);
}

void CodonM2aModel::CollectLengthSuffStat()    {

    lengthsuffstatarray->Clear();
    phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
}

void CodonM2aModel::MoveLambda()	{

    lambdasuffstat.Clear();
    branchlength->AddSuffStat(lambdasuffstat);
    MoveLambda(1.0,10);
    MoveLambda(0.3,10);
    blhypermean->SetAllBranches(1.0/lambda);
}

double CodonM2aModel::MoveLambda(double tuning, int nrep)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - LambdaHyperLogPrior() - LambdaHyperSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        lambda *= e;
        deltalogprob += LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
        deltalogprob += m;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            lambda /= e;
        }
        ntot++;
    }
    return nacc/ntot;
}

//
// Omega mixture 
//

void CodonM2aModel::CollectPathSuffStat()	{

    sitepathsuffstatarray->Clear();
    phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
}

void CodonM2aModel::CollectComponentPathSuffStat()	{

    componentpathsuffstatarray->Clear();
    sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
}

void CodonM2aModel::MoveOmega() 	{

    CollectOmegaSuffStat();

    OmegaHyperSlidingMove(purom,0.1,10,0,1);
    OmegaHyperSlidingMove(purw,1,10,0,1);
    if (pi != 0)    {
        OmegaHyperScalingMove(dposom,1,10);
        OmegaHyperSlidingMove(posw,1,10,0,1);
    }
    if ((pi != 0) && (pi != 1))    {
        SwitchPosWeight(10);
    }

    ResampleAlloc();
}

void CodonM2aModel::CollectOmegaSuffStat()	{

    siteomegasuffstatarray->Clear();
    siteomegasuffstatarray->AddSuffStat(*sitecodonmatrixarray,*sitepathsuffstatarray);
}

void CodonM2aModel::ResampleAlloc()	{
    OmegaSuffStatLogProb();
    sitealloc->GibbsResample(sitepostprobarray);
}

double CodonM2aModel::OmegaHyperScalingMove(double& x, double tuning, int nrep)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - OmegaLogPrior() - OmegaSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        x *= e;
        deltalogprob += OmegaLogPrior() + OmegaSuffStatLogProb();
        deltalogprob += m;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            x /= e;
        }
        ntot++;
    }
    return nacc/ntot;
}

double CodonM2aModel::OmegaHyperSlidingMove(double& x, double tuning, int nrep, double min, double max)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - OmegaLogPrior() - OmegaSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        x += m;
        if (max > min)  {
            while ((x < min) || (x > max))  {
                if (x < min)    {
                    x = 2*min - x;
                }
                if (x > max)    {
                    x = 2*max - x;
                }
            }
        }
        deltalogprob += OmegaLogPrior() + OmegaSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            x -= m;
        }
        ntot++;
    }
    return nacc/ntot;
}

double CodonM2aModel::DrawBetaPosWeight()    {
    double alpha = poswhypermean / poswhyperinvconc;
    double beta = (1-poswhypermean) / poswhyperinvconc;
    double a = Random::sGamma(alpha);
    double b = Random::sGamma(beta);
    double ret = a / (a+b);
    return ret;
}

double CodonM2aModel::SwitchPosWeight(int nrep)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double bkposw = posw;
        double deltalogprob = - PosSwitchLogProb() - OmegaSuffStatLogProb();
        if (posw)   {
            posw = 0;
        }
        else    {
            posw = DrawBetaPosWeight();
        }
        deltalogprob += PosSwitchLogProb() + OmegaSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            posw = bkposw;
        }
        ntot++;
    }
    return nacc/ntot;
}

//
// nucleotide parameters
//

void CodonM2aModel::CollectNucPathSuffStat()   {
    UpdateMatrices();
    nucpathsuffstat.Clear();
    nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
}

void CodonM2aModel::MoveNucRates()	{

    CollectComponentPathSuffStat();
    CollectNucPathSuffStat();

    MoveRR(0.1,1,3);
    MoveRR(0.03,3,3);
    MoveRR(0.01,3,3);

    MoveNucStat(0.1,1,3);
    MoveNucStat(0.01,1,3);

    UpdateMatrices();
}

double CodonM2aModel::MoveRR(double tuning, int n, int nrep)	{
    double nacc = 0;
    double ntot = 0;
    double bk[Nrr];
    for (int rep=0; rep<nrep; rep++)	{
        for (int l=0; l<Nrr; l++)	{
            bk[l] = nucrelrate[l];
        }
        double deltalogprob = -NucRatesLogPrior() - NucRatesSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(nucrelrate,Nrr,tuning,n);
        deltalogprob += loghastings;
        UpdateNucMatrix();
        deltalogprob += NucRatesLogPrior() + NucRatesSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            for (int l=0; l<Nrr; l++)	{
                nucrelrate[l] = bk[l];
            }
            UpdateNucMatrix();
        }
        ntot++;
    }
    return nacc/ntot;
}

double CodonM2aModel::MoveNucStat(double tuning, int n, int nrep)	{
    double nacc = 0;
    double ntot = 0;
    double bk[Nnuc];
    for (int rep=0; rep<nrep; rep++)	{
        for (int l=0; l<Nnuc; l++)	{
            bk[l] = nucstat[l];
        }
        double deltalogprob = -NucRatesLogPrior() - NucRatesSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(nucstat,Nnuc,tuning,n);
        deltalogprob += loghastings;
        UpdateNucMatrix();
        deltalogprob += NucRatesLogPrior() + NucRatesSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            for (int l=0; l<Nnuc; l++)	{
                nucstat[l] = bk[l];
            }
            UpdateNucMatrix();
        }
        ntot++;
    }
    return nacc/ntot;
}


// summary statistics

double CodonM2aModel::GetTotalLength()	{
    double tot = 0;
    for (int j=1; j<Nbranch; j++)	{
        tot += branchlength->GetVal(j);
    }
    return tot;
}

double CodonM2aModel::GetMeanOmega()   {
    return posw*(1 + dposom) + (1-posw)*(purw*purom + (1-purw));
}

double CodonM2aModel::GetEntropy(const std::vector<double>& profile, int dim) const {
    double tot = 0;
    for (int i=0; i<dim; i++)	{
        tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
    }
    return tot;
}

void CodonM2aModel::TraceHeader(std::ostream& os)  {
    os << "#logprior\tlnL\tlength\t";
    os << "purom\tposom\tpurw\tposw\t";
    os << "statent\t";
    os << "rrent\n";
}

void CodonM2aModel::Trace(ostream& os) {	
    os << GetLogPrior() << '\t';
    os << GetLogLikelihood() << '\t';
    os << GetTotalLength() << '\t';
    os << purom << '\t' << dposom+1 << '\t' << purw << '\t' << posw << '\t';
    os << GetEntropy(nucstat,Nnuc) << '\t';
    os << GetEntropy(nucrelrate,Nrr) << '\n';
    SubMatrix::diagerr = 0;
}

void CodonM2aModel::TracePostProb(ostream& os) {
    for (int i=0; i<GetNsite(); i++)    {
        os << sitepostprobarray[i][2] << '\t';
    }
    os << '\n';
}

