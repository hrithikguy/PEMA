#if !defined(CLUSTER_H)
#define  CLUSTER_H

#define CRI(a) ((a)==0?"BIC":(a)==1?"AIC":"AICc")

#define NAME "MACML"	//Program name
#define FULLNAME "Model Averaging Clustering by Maximum Likelihood"
#define FUNCTION "Cluster sequences into heterogeneous regions with specified site types"
#define CITATION "Zhang, Z., and Townsend, J.P. (2009) Maximum-likelihood model averaging to profile clustering of site types across discrete linear sequences. PLoS Computational Biology 5(6): e1000421."	

#define VERSION "Version: 1.1.2 [December 8, 2011]"
/**************** Updates for 1.1.1: **********************
   * Optimize the code
	* Add a jump table to reduce the number of clustering
	* Add more options for input parameters
*********************************************************/

//#define VERSION "Version: 1.0 (January 07, 2009)"
/*** Initial versioin ***/

#include "base.h"
#include <map>

using namespace std;

class CandidateModels {

public:
	double CW; //criterion and weight
	
	long x_pos_start, x_pos_end, y_pos_start, y_pos_end;
	double p0, pc;
	long x_cs, x_ce, y_cs, y_ce;

	double InL0, InL;
	double AIC0, AIC;
	double AICc0, AICc;
	double BIC0, BIC;

	CandidateModels() {}

	//For all possible candidate models
	CandidateModels(double CW, long x_pos_start, long x_pos_end, long y_pos_start, long y_pos_end, long x_cs, long x_ce, long y_cs, long y_ce, double p0, double pc){
		this->CW = CW;
		this->x_pos_start = x_pos_start;		
		this->x_pos_end = x_pos_end;
		this->x_cs = x_cs;		
		this->x_ce = x_ce;
		this->y_pos_start = y_pos_start;		
		this->y_pos_end = y_pos_end;
		this->y_cs = y_cs;		
		this->y_ce = y_ce;
		this->p0 = p0;		
		this->pc = pc;
	}

	//For best-fit selected models
	CandidateModels(long x_pos_start, long x_pos_end, long y_pos_start, long y_pos_end, long x_cs, long x_ce, long y_cs, long y_ce, double p0, double pc, double InL0, double InL,
		double AIC0, double AIC, double AICc0, double AICc, double BIC0, double BIC) {		
		this->x_pos_start = x_pos_start;		
		this->x_pos_end = x_pos_end;
		this->x_cs = x_cs;		
		this->x_ce = x_ce;
		this->y_pos_start = y_pos_start;		
		this->y_pos_end = y_pos_end;
		this->y_cs = y_cs;		
		this->y_ce = y_ce;
		this->p0 = p0;		
		this->pc = pc;
		this->InL0 = InL0;		
		this->InL = InL;
		this->AIC0 = AIC0;		
		this->AIC = AIC;
		this->AICc0 = AICc0;		
		this->AICc = AICc;
		this->BIC0 = BIC0;		
		this->BIC = BIC;
	}

};

class CI {//confidence interval

public:
	double weight; //weight
	double p;	//probability

	CI() {}

	CI(double weight, double p) {
		this->weight = weight;	
		this->p = p;
	}
};

class Cluster: public Base {
	
public:	
	Cluster();
	~Cluster();	

	int Run(int argc, const char*argv[]);
	//Clustering using maximum likelihood

	int RunML(vector <double> patient_info, vector < vector <double> > seq2d);

	//int RunML(vector<vector <double> > seq2d);
	
	//Cluster (sub)sequences based on divide-and-conquer

	int ClusterSubSeq(vector <double> patient_info, vector <vector <double> > seq2d, long x_pos_start, long x_pos_end, long y_pos_start, long y_pos_end);

	//int ClusterSubSeq(vector <vector <double> > seq2d, long x_pos_start, long x_pos_end, long y_pos_start, long y_pos_end);

	//int ClusterSubSeq(string seq, long pos_start, long pos_end);

	//Model averaging based on all possible candidate models
	int ModelAveraging(long x_pos_start, long x_pos_end, long y_pos_start, long y_pos_end, long x_cs_max, long x_ce_max, long y_cs_max, long y_ce_max, double p0_max, double pc_max, double MIN_cri);

	//int ModelAveraging(long pos_start, long pos_end, long cs, long ce, double p0, double pc, double MIN_cri);

	//Initialize
	int init(vector<vector <double> > seq2d);

	//Output results
	int output(vector<vector <double> > seq2d);

protected:
	double BinomialProb(long n, long i, double p);
	//double BinomialProb(double n, double i, double p);
	double getp0pc(vector <vector <double> > seq2d, long x_pos_start, long x_pos_end, long x_cs, long x_ce, long y_pos_start, long y_pos_end, long y_cs, long y_ce, double &p0, double &pc, long &ni, long &no, long &in_total, long &out_total);

	//double getp0pc(string seq, long pos_start, long pos_end, long cs, long ce, double &p0, double &pc, long &ns, long &nc, long &ne);
	//Jump table for cs, ce, to reduce the number of clustering
	int jumpTable(string seq, long start_pos, long end_pos);
	//Parse parameters
	bool parseParameter(int argc, const char* argv[]);
	
public:
	//A vector of selected models
	vector<CandidateModels> vec_SelectedModels;

	  //Count the number of clusters
	  int flag_found_cluster;


	//H=Hot Spot, C=Cold Spot, default=-
	vector<vector <char> > vec_spot;
	//rate by model selection
	vector<vector <double> > vec_MS_rate;
	//rate by model averaging
	vector<vector <double> >vec_MA_rate;
	//Lower CI
	vector<vector <double> > vec_lower_rate;
	//Upper CI
	vector<vector <double> > vec_upper_rate;

protected:
	//Name of input sequence file
	string seqfile;
	//vector of sequence
	vector<string> seq;

	vector <double> patient_info;
	vector <vector <double> > seq2d;

	//vector of sequence name
	vector<string> seqname;

	//Confidence interval for model averaging, default=0.95
	double confidence_interval;

private:
	//Sequence type, default=0
	int seq_type;
	//Criterion, 0=BIC, 1=AIC, 2=AICc, default=0
	int criterion_type;
	int sample_size;
	//Only model selection, default=0
	int MS_only;

	//Jump table for ce
	map <long, long> CETABLE;
	//Jump table for cs
	map <long, long> CSTABLE;

	//A vector of all candidate models
	  //For short definition
	vector<CandidateModels> vec_AllModels;
	//For the long definition
	vector<CandidateModels> vec_AllModelsFull;
};


#endif


