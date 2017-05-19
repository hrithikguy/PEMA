#include "HC.h"


Cluster::Cluster() {
	patient_info.clear();
	seq2d.clear();

	seqfile = "";

	//seq_type = 0;
	criterion_type = 0;
	MS_only = 0;
	sample_size = 500;
	  
	confidence_interval = 0.95;
	confidence_interval = (1.0 - confidence_interval)/2.0;

	flag_found_cluster=0;
}

Cluster::~Cluster() {
	patient_info.clear();
	seq2d.clear();

	vec_lower_rate.clear();
	vec_upper_rate.clear();
	vec_MS_rate.clear();
	vec_MA_rate.clear();
	vec_spot.clear();

	vec_SelectedModels.clear();
	vec_AllModels.clear();
	  vec_AllModelsFull.clear();
}


double Cluster::BinomialProb(long n, long i, double p) {
	double prob = (p<SMALLVALUE || i==0)?0:i*log(p);
	prob += (1.0-p<SMALLVALUE || n-i==0)?0:(n-i)*log(1.0-p);
	return prob;
}

double Cluster::getp0pc(vector <vector <double> > seq2d, long x_pos_start, long x_pos_end, long x_cs, long x_ce, long y_pos_start, long y_pos_end, long y_cs, long y_ce, double &p0, double &pc, long &ni, long &no, long &in_total, long &out_total) {
	
	// long j;	

	// ns=nc=ne=0;	
	// for (j=pos_start; j<=cs-1; j++) if (seq[j]=='1') ns++;
	// for (j=cs;  j<=ce;   j++) if (seq[j]=='1') nc++;
	// for (j=ce+1;j<=pos_end;  j++) if (seq[j]=='1') ne++;
	
	// pc = (double)nc/(ce-cs+1);
	// if (cs==pos_start && pos_end==ce) {
	// 	p0 = 0.0;
	// }
	// else {
	// 	p0 = (double)(ns+ne)/(cs-pos_start+pos_end-ce);
	// }

	long j,k;
	ni = no = 0;

	in_total = 0;
	out_total = 0;

	for (j = x_cs; j <= x_ce; ++j) {
		for (k = y_cs; k <= y_ce; ++k) {
			if (seq2d[k][j] == 1) {
				//cout << "FOUND SOMETHING INSIDE A CLUSTER!!!!!!!!!!!!!!!!!!!!!" << endl;
				++ni;
				++in_total;
			}
			if (seq2d[k][j] == 0) {
				++in_total;
			}
		}
	}

	if (in_total < 2) {
		return 0;
	}


	for (j = x_pos_start; j <= x_pos_end; ++j) {
		for (k = y_pos_start; k <= y_pos_end; ++k) {
			if (seq2d[k][j] == 1) {
				//cout << "FOUND SOMETHING INSIDE A CLUSTER!!!!!!!!!!!!!!!!!!!!!" << endl;
				++no;
				++out_total;
			}
			if (seq2d[k][j] == 0) {
				++out_total;
			}
		}
	}


	out_total = out_total - in_total;
	no = no - ni;


	if (in_total == 0) {
		pc = 0.0;
	}
	else {
		pc = (double) ni/in_total;
	}

	//if (x_pos_start == 0 && x_pos_end == seq2d[0].size()-1 && y_pos_start == 0 && y_pos_end == seq2d.size() - 1) {
	if (out_total == 0) {
		p0 = 0.0;
		cout << "setting p0 to 0 automatically" << endl;
	}
	else {
		p0 = (double) no/ (double) out_total;
	}


	return 1;
}



// int Cluster::jumpTable(vector <double> patient_info, long start_pos, long end_pos) {

// 	long i, j;

// 	CSTABLE.clear();
// 	CETABLE.clear();
	
// 	for (i=start_pos; i<=end_pos; i++) {		
// 		if(MS_only==0) {
// 			CSTABLE[i] = CETABLE[i] = i+1;
// 		}
// 		else {
// 			j = i + 1;
// 			while (seq[i]==seq[j] && j<=end_pos) {
// 				j++;
// 			}		
// 			if (seq[i]!=seq[j] && j!=i+1) {
// 				CSTABLE[i] = j;
// 				CETABLE[i] = j-1;
// 			}
// 			else {
// 				CSTABLE[i] = i+1;
// 				CETABLE[i] = i+1;
// 			}
// 		}
// 	}

// 	return 1;
// }

int Cluster::Run(int argc, const char*argv[]) {
	
	try {
		cout<<NAME<<" ("<<FULLNAME<<"), "<<VERSION<<endl;

		if (parseParameter(argc, argv)==false) {			
			cout<<"Function: "<<FUNCTION<<endl;
			cout<<"Usage: "<<NAME<<" [arguments]"<<endl<<endl;
			cout<<"Arguments:"<<endl;
			cout<<"  -i\tInput sequence file name [string, required]"<<endl;			
			cout<<"  -c\tCriterion for clustering [integer, optional], {0:BIC || 1:AIC || 2:AICc}, default = 0"<<endl;
			cout<<"  -t\tType of input sequence [integer, optional], {0: aligned sequences || 1: a vector of 0's and 1's], default = 0"<<endl;
			cout<<"  -m\tModel selection and model averaging  [integer, optional], {0: use both model selection and model averaging || 1: use only model selection], default = 0"<<endl;
			cout << "-s Number of samples of models during model averaging. Default is 500." << endl;
			cout<<endl<<"Please send suggestions or bugs to Zhang Zhang at zhang.zhang@yale.edu."<<endl<<endl;
			throw 1;
		}
		
		// if (readFasta(seqfile, seqname, seq)==0) throw "reading fasta sequences";

		read_input(seqfile, patient_info, seq2d);
		// int jj, kk;
		// for (jj = 0; jj < patient_info.size(); ++jj) {
		// 	cout << patient_info[jj] << endl;
		// 	for (kk = 0; kk < seq2d[jj].size(); ++kk) {
		// 		cout << seq2d[jj][kk] << endl;
		// 	}
		// }
		
		cout<<endl<<"Reference: "<<CITATION<<endl;

		cout<<endl<<"Criterion used: "<<CRI(criterion_type)<<endl<<endl;

		static time_t time_start = time(NULL);
		
		// int i;
		// if (seq_type==0) {//a set of aligned sequences
		// 	for (i=0; i<seq.size(); i++) {
		// 		cout<<seqname[i].c_str()<<"\t"<<seq[i].c_str()<<endl;
		// 		if (seq[0].length()!=seq[i].length()) throw "unequal length of input sequences";
		// 	}
		// 	string consensus = generalInput(seq);
		// 	cout<<"Consensus"<<"\t"<<consensus.c_str()<<endl;
		// 	RunML(consensus);			
		// 	cout<<endl;
		// }
		// else if (seq_type==1) {//a vector of 0's and 1's
		// 	for (i=0; i<seq.size(); i++) {				
		// 		cout<<">"<<seqname[i].c_str()<<endl;
		// 		RunML(seq[i]);				
		// 	}
		// }

		RunML(patient_info, seq2d);

		//display on screen
		cout<<"Mission accomplished. (Time elapsed: ";
		time_t t = time(NULL)-time_start;
		int h=t/3600, m=(t%3600)/60, s=t-(t/60)*60;
		if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
		else   cout<<m<<":"<<s<<")"<<endl;		
	}
	catch (int exp) {}
	catch (const char *e) {
		cout<<"Error: "<<e<<endl;
	}

	return 1;
}

int Cluster::init(vector <vector <double> > seq2d) {
	long X_len = seq2d[0].size();
	long Y_len = seq2d.size();
	

	long i;
	vector <double> zeros(X_len, 0.0);
	vector <char> dashes(X_len, '-');

	for (i = 0; i < Y_len; ++i) {
		vec_MS_rate.push_back(zeros);
		vec_MA_rate.push_back(zeros);
		vec_lower_rate.push_back(zeros);
		vec_upper_rate.push_back(zeros);
		vec_spot.push_back(dashes);
	}
	// vec_MS_rate.resize(N,0.0);
	// vec_spot.resize(N,'-');
	// vec_MA_rate.resize(N,0.0);
	// vec_lower_rate.resize(N,0.0);
	// vec_upper_rate.resize(N,0.0);
	vec_SelectedModels.clear();
	vec_AllModels.clear();
	//jumpTable(seq, 0, N-1);

	 vec_AllModelsFull.clear();
	return 1;
}

int Cluster::output(vector< vector <double> >  seq2d) {

	long i,j;
	  cout<<endl<<"//Results based on model selection:"<<endl;
	  if(vec_SelectedModels.size()==0){
	    cout<<"Note: No clusters are found in the sequence!"<<endl;
	    cout<<"There are not enough variants!"<<endl;
	  }

	for(i=0; i<vec_SelectedModels.size(); i++) {
	  if(flag_found_cluster>0){
	    cout<<vec_SelectedModels[i].x_pos_start<<" ~ "<<vec_SelectedModels[i].x_pos_end << ", ";
	   	cout<<vec_SelectedModels[i].y_pos_start<<" ~ "<<vec_SelectedModels[i].y_pos_end;

	    cout<<"\tx_cs= "<<vec_SelectedModels[i].x_cs<<"\tx_ce= "<<vec_SelectedModels[i].x_ce << ", ";
	    cout<<"\ty_cs= "<<vec_SelectedModels[i].y_cs<<"\ty_ce= "<<vec_SelectedModels[i].y_ce;

	    cout<<"\tp0= "<<vec_SelectedModels[i].p0<<"\tpc= "<<vec_SelectedModels[i].pc;
	    cout<<"\tInL0= "<<vec_SelectedModels[i].InL0<<"\tInL= "<<vec_SelectedModels[i].InL;
	    cout<<"\tAIC0= "<<vec_SelectedModels[i].AIC0<<"\tAIC= "<<vec_SelectedModels[i].AIC;
	    cout<<"\tAICc0= "<<vec_SelectedModels[i].AICc0<<"\tAICc= "<<vec_SelectedModels[i].AICc;
	    cout<<"\tBIC0= "<<vec_SelectedModels[i].BIC0<<"\tBIC= "<<vec_SelectedModels[i].BIC;
	    cout<<endl;
	  }else{
	    if(i==0){
	      cout<<"Note: No clusters are found in the sequence!"<<endl;
	      cout<<"Null Model result:"<<endl;
	      cout<<vec_SelectedModels[i].x_pos_start<<" ~ "<<vec_SelectedModels[i].x_pos_end << ", ";
	   	cout<<vec_SelectedModels[i].y_pos_start<<" ~ "<<vec_SelectedModels[i].y_pos_end;

	    cout<<"\tx_cs= "<<vec_SelectedModels[i].x_cs<<"\tx_ce= "<<vec_SelectedModels[i].x_ce << ", ";
	    cout<<"\ty_cs= "<<vec_SelectedModels[i].y_cs<<"\ty_ce= "<<vec_SelectedModels[i].y_ce;

	    cout<<"\tp0= "<<vec_SelectedModels[i].p0<<"\tpc= "<<vec_SelectedModels[i].pc;
	    cout<<"\tInL0= "<<vec_SelectedModels[i].InL0<<"\tInL= "<<vec_SelectedModels[i].InL;
	    cout<<"\tAIC0= "<<vec_SelectedModels[i].AIC0<<"\tAIC= "<<vec_SelectedModels[i].AIC;
	    cout<<"\tAICc0= "<<vec_SelectedModels[i].AICc0<<"\tAICc= "<<vec_SelectedModels[i].AICc;
	    cout<<"\tBIC0= "<<vec_SelectedModels[i].BIC0<<"\tBIC= "<<vec_SelectedModels[i].BIC;
	      cout<<endl<<endl;
	    }
	    if(i==1){
	      cout<<"The cluster with the best model for the sequence (not significant):"<<endl;
	    

              cout<<vec_SelectedModels[i].x_pos_start<<" ~ "<<vec_SelectedModels[i].x_pos_end << ", ";
	   	cout<<vec_SelectedModels[i].y_pos_start<<" ~ "<<vec_SelectedModels[i].y_pos_end;

	    cout<<"\tx_cs= "<<vec_SelectedModels[i].x_cs<<"\tx_ce= "<<vec_SelectedModels[i].x_ce << ", ";
	    cout<<"\ty_cs= "<<vec_SelectedModels[i].y_cs<<"\ty_ce= "<<vec_SelectedModels[i].y_ce;

	    cout<<"\tp0= "<<vec_SelectedModels[i].p0<<"\tpc= "<<vec_SelectedModels[i].pc;
	    cout<<"\tInL0= "<<vec_SelectedModels[i].InL0<<"\tInL= "<<vec_SelectedModels[i].InL;
	    cout<<"\tAIC0= "<<vec_SelectedModels[i].AIC0<<"\tAIC= "<<vec_SelectedModels[i].AIC;
	    cout<<"\tAICc0= "<<vec_SelectedModels[i].AICc0<<"\tAICc= "<<vec_SelectedModels[i].AICc;
	    cout<<"\tBIC0= "<<vec_SelectedModels[i].BIC0<<"\tBIC= "<<vec_SelectedModels[i].BIC;


              cout<<endl;
	    }
	  }
	}
	
	if (MS_only==0) {
		cout<<endl<<"//Results based on model averaging:"<<endl;
		cout.setf(ios::left); 
		int width=15;
		cout.width(width); cout<<"X Position";
		cout.width(width); cout<<"Y Position";
		cout.width(width); cout<<"Site";
		cout.width(width); cout<<"MS";
		cout.width(width); cout<<"MA";
		cout.width(width); cout<<"Lower_CI";
		cout.width(width); cout<<"Upper_CI";
		cout.width(width); cout<<"Hot|Cold";
		cout<<endl;
		
		for(i=0; i<seq2d.size(); i++) {	
			for (j = 0; j < seq2d[0].size(); j++) {
				cout.width(width);cout<<j<< i << "\t";
				cout.width(width);cout<<seq2d[i][j]<<"\t";
				cout.width(width);cout<<vec_MS_rate[i][j]<<"\t";
				cout.width(width);cout<<vec_MA_rate[i][j]<<"\t";
				cout.width(width);cout<<vec_lower_rate[i][j]<<"\t";
				cout.width(width);cout<<vec_upper_rate[i][j]<<"\t";
				cout.width(width);cout<<vec_spot[i][j];
				cout<<endl;
			}
		}

		cout<<endl<<"Abbreviation:  MS=Model Selection; MA=Model Averaging; CI=Confidence Interval; H=Hot; C=Cold;"<<endl;
	}
	
	cout<<endl<<"#End of clustering"<<endl<<endl;

	return 1;
}


int Cluster::RunML(vector <double> patient_info, vector < vector <double> > seq2d) {

	cout << "Starting RunML " << endl;
	//Initialize
	init(seq2d);
	cout << "Initialized seq2d" << endl;

	//Cluster sequences based on maximum likelihood
	ClusterSubSeq(patient_info, seq2d, 0, seq2d[0].size()-1, 0, seq2d.size()-1);

	cout << "Finished clustering " << endl;

	//output results
	output(seq2d);

	return 1;
}


int Cluster::ClusterSubSeq(vector <double> patient_info, vector <vector <double> > seq2d, long x_pos_start, long x_pos_end, long y_pos_start, long y_pos_end) {
	

	// long test1, test2;
	// for (test1 = 0; test1 < seq2d.size(); ++test1) {
	// 	for (test2 =0 ; test2 < seq2d[test1].size(); ++test2) {
	// 		cout << seq2d[test1][test2] << "\t";
	// 	}
	// 	cout << endl;
	// }

	long X_size = x_pos_end - x_pos_start + 1;
	long Y_size = y_pos_end - y_pos_start + 1;
	long i,j, n=0;
	for (i=x_pos_start; i<=x_pos_end; i++) {
		for (j = y_pos_start; j <= y_pos_end; j++) {

			if (seq2d[j][i]== 1) n++;
		}
	}

	//cout << "bp1" << endl;
	
	if (X_size <= 2 || Y_size <= 2 || n==0 || n== X_size * Y_size) {
		return 1;
	}
	
	double AIC0, AICc0, BIC0, InL0 = BinomialProb(X_size * Y_size, n, double(n)/(X_size * Y_size));
	AIC0 = AICc0 = BIC0 = -2*InL0; //parameter = 0 for null model

	double InL = InL0;
	double AIC = AIC0;
	double AICc = AICc0;
	double BIC = BIC0;
	
	long x_cs, x_ce, y_cs, y_ce;
	double p0_max, pc_max, x_cs_max, x_ce_max, y_cs_max, y_ce_max;
	int isFound=0; //Whether found the lowest AIC/BIC		
	
	double MIN_cri = 1000000;
	double para = 4.0;

	  vec_AllModels.clear();
	  vec_AllModelsFull.clear();

	//cout << "bp2" << endl;


	for (x_cs = x_pos_start; x_cs <= x_pos_end; ++x_cs) {
			//cout << "bp3" << endl;


		for (y_cs = y_pos_start; y_cs <= y_pos_end; ++y_cs) {
			if (y_cs > 0 && patient_info[y_cs-1] == patient_info[y_cs]) {
				continue;
			}
			for(x_ce=x_cs+1; x_ce <= x_pos_end; ++x_ce) {
				for (y_ce = y_cs + 1; y_ce <= y_pos_end; ++y_ce) {
					if (y_ce < seq2d.size() - 1 && patient_info[y_ce] == patient_info[y_ce + 1]) {
						continue;
					}
					//cout << "trying a cluster" << endl;
					//cout << x_cs << "\t" << x_ce << "\t" << y_cs << "\t" << y_ce << endl;
					//cout << "that was a cluster" << endl;
					//if ((cs==pos_start && ce==pos_end) || seq[cs]!=seq[ce]) continue;			
					if (x_cs==x_pos_start && x_ce==x_pos_end && y_cs==y_pos_start && y_ce==y_pos_end ) {
						para = 0.0;			
					}
					else {
						if (x_pos_end-x_ce<=1 || x_cs-x_pos_start<=1 || y_pos_end-y_ce<=1 || y_cs-y_pos_start<=1) continue;
						para = 4.0;
					}

					long ni, no;
					long in_total, out_total;
					double p0, pc;
					double at_least_two_data_points;
					at_least_two_data_points = getp0pc(seq2d, x_pos_start, x_pos_end, x_cs, x_ce, y_pos_start, y_pos_end, y_cs, y_ce, p0, pc, ni, no, in_total, out_total);
					
					if (at_least_two_data_points == 0) {
						continue;
					}
					//cout << p0 << "\t" << pc << endl;
					double InL_tmp = BinomialProb(in_total, ni, pc) + BinomialProb(out_total, no, p0);
					
					//cout << InL_tmp << endl;
					double AIC_tmp  = -2*InL_tmp + 2*para;			
					double AICc_tmp = AIC_tmp;
					if (X_size * Y_size-para-1>0.0) AICc_tmp += 2*para*(para+1)/(X_size * Y_size-para-1);
					else AICc_tmp = 2*AIC_tmp;//sequence too short, give more penalty
					double BIC_tmp = -2*InL_tmp + para*log(double(X_size * Y_size));
					
					double cri, cri0;
					if (criterion_type==0) {//BIC
						cri0 = BIC;
						cri = BIC_tmp;
					}
					else if (criterion_type==1) {//AIC
						cri0 = AIC;
						cri = AIC_tmp;
					}
					else if (criterion_type==2) {//AICc
						cri0 = AICc;
						cri = AICc_tmp;
					}

					  //cout<<"cri="<<cri<<";"<<pos_start<<";"<<pos_end<<";"<<cs<<";"<<ce<<";"<<p0<<";"<<pc<<endl;
					  if (MIN_cri > cri){
					    MIN_cri = cri;
					  }




					CandidateModels tmp_L(cri, x_pos_start, x_pos_end, y_pos_start, y_pos_end, x_cs, x_ce, y_cs, y_ce, p0, pc);
					vec_AllModels.push_back(tmp_L);
					//cout << "adding to vec_models" << endl;
					  //If there is no cluster, we use this vector to find out the best model except the null model.
					  CandidateModels tmp_L_full(x_pos_start, x_pos_end, y_pos_start, y_pos_end, x_cs, x_ce, y_cs, y_ce, p0, pc, InL0, InL_tmp, AIC0, AIC_tmp, AICc0, AICc_tmp, BIC0, BIC_tmp);
					  vec_AllModelsFull.push_back(tmp_L_full);

					if (cri<cri0) {
						if (x_cs-x_pos_start>1 && x_pos_end-x_ce>1 && y_cs-y_pos_start>1 && y_pos_end-y_ce>1) {
							//TODO: ask about this
							//if ( (p0<pc && seq[cs]=='1' && seq[cs-1]=='0' && seq[ce+1]=='0') || (p0>pc && seq[cs]=='0' && seq[cs-1]=='1' && seq[ce+1]=='1') ) {
							isFound = 1;
							
							x_cs_max = x_cs;
							x_ce_max = x_ce;
							y_cs_max = y_cs;
							y_ce_max = y_ce;

							p0_max = p0;
							pc_max = pc;
							
							InL = InL_tmp;
							AIC = AIC_tmp;
							AICc = AICc_tmp;
							BIC = BIC_tmp;


							cout<<tmp_L.x_pos_start<<" ~ "<<tmp_L.x_pos_end << ", ";
						   	cout<<tmp_L.y_pos_start<<" ~ "<<tmp_L.y_pos_end;

						    cout<<"\tx_cs= "<<tmp_L.x_cs<<"\tx_ce= "<<tmp_L.x_ce << ", ";
						    cout<<"\ty_cs= "<<tmp_L.y_cs<<"\ty_ce= "<<tmp_L.y_ce;

						    cout<<"\tp0= "<<tmp_L.p0<<"\tpc= "<<tmp_L.pc;
						    cout<<"\tInL0= "<<tmp_L.InL0<<"\tInL= "<<tmp_L.InL;
						    cout<<"\tAIC0= "<<tmp_L.AIC0<<"\tAIC= "<<tmp_L.AIC;
						    cout<<"\tAICc0= "<<tmp_L.AICc0<<"\tAICc= "<<tmp_L.AICc;
						    cout<<"\tBIC0= "<<tmp_L.BIC0<<"\tBIC= "<<tmp_L.BIC;
						    cout<<endl;

						    long test1, test2;
							for (test2 = y_cs_max ; test2 <= y_ce_max; ++test2) {
								for (test1 = x_cs_max; test1 <= x_ce_max; ++test1) {

									cout << seq2d[test2][test1];
								}
								cout << endl;
							}
						    

							//}
						}
		 			}
		 		}
			}
		}
	}
		

			//cout << "bp4" << endl;
	
	  if (isFound==0) {
	    //If there is no cluster, we keep the null model, and the best model which is not significant.
	    if(flag_found_cluster==0){
	      //Keep the null model
	      double p_tmp=(double)n/(double)(X_size*Y_size);
	      CandidateModels nullmodel(x_pos_start, x_pos_end, y_pos_start, y_pos_end,  x_pos_start, x_pos_end, y_pos_start, y_pos_end, p_tmp , p_tmp, InL0, InL, AIC0, AIC, AICc0, AICc, BIC0, BIC);
	      vec_SelectedModels.push_back(nullmodel);

	      	//cout << "bp4.1" << endl;

	      //Find out the smallest cri except InL0
	      //vec_AllModelsFull
	      long smallest_j=0;
	      double smallest_cri=100000;
	      for(long j=0;j<vec_AllModelsFull.size();j++){
			if(criterion_type==0 && vec_AllModelsFull[j].BIC < smallest_cri && vec_AllModelsFull[j].BIC != BIC){
			  smallest_cri=vec_AllModelsFull[j].BIC;
			  smallest_j=j;
			}else if (criterion_type==1 && vec_AllModelsFull[j].AIC < smallest_cri && vec_AllModelsFull[j].AIC != AIC){
			  smallest_cri=vec_AllModelsFull[j].AIC;
			  smallest_j=j;
			}else if (criterion_type==2 && vec_AllModelsFull[j].AICc < smallest_cri && vec_AllModelsFull[j].AICc != AICc){
	                  smallest_cri=vec_AllModelsFull[j].AICc;
	                  smallest_j=j;
			} 
	      }
	      	      	//cout << "bp4.2" << endl;

	      	      	cout << vec_AllModelsFull.size() << endl;
	      CandidateModels bestmodel=vec_AllModelsFull[smallest_j];
	      	      	      	//cout << "bp4.21" << endl;

	      vec_SelectedModels.push_back(bestmodel);
	      	      	      	//cout << "bp4.3" << endl;

	    }
	    	      	      	//cout << "bp4.4" << endl;

	    return 1;
	  }else{
	    flag_found_cluster++;
	  }

	  	//cout << "bp5" << endl;


	CandidateModels selectedmodel(x_pos_start,x_pos_end,y_pos_start, y_pos_end, x_cs_max,x_ce_max,y_cs_max, y_ce_max,p0_max,pc_max,InL0,InL,AIC0,AIC,AICc0,AICc,BIC0,BIC);
	vec_SelectedModels.push_back(selectedmodel);

	if (MS_only==0) ModelAveraging(x_pos_start, x_pos_end, y_pos_start, y_pos_end, x_cs_max, x_ce_max, y_cs_max, y_ce_max, p0_max, pc_max, MIN_cri);

	/* Divide and Conquer */
	if (x_ce_max!=x_pos_end || x_cs_max!=x_pos_start || y_ce_max!=y_pos_end || y_cs_max!=y_pos_start ) {				
		//if (cs_max>pos_start+1) ClusterSubSeq(seq, pos_start, cs_max-1);
		//if (x_ce_max<pos_end-1) ClusterSubSeq(seq, ce_max+1, pos_end);
		if (x_cs_max > x_pos_start +1 && y_cs_max > y_pos_start + 1) ClusterSubSeq(patient_info, seq2d, x_pos_start, x_cs_max, y_pos_start, y_cs_max);
		if (y_cs_max > y_pos_start + 1) ClusterSubSeq(patient_info, seq2d, x_cs_max, x_ce_max, y_pos_start, y_cs_max);
		if (x_ce_max < x_pos_end - 1 && y_cs_max > y_pos_start + 1) ClusterSubSeq(patient_info, seq2d, x_ce_max, x_pos_end, y_pos_start, y_cs_max);
		if (x_cs_max > x_pos_start + 1) ClusterSubSeq(patient_info, seq2d, x_pos_start, x_cs_max, y_cs_max, y_ce_max);
		if (x_cs_max<x_ce_max && y_cs_max < y_ce_max) ClusterSubSeq(patient_info, seq2d, x_cs_max, x_ce_max, y_cs_max, y_ce_max);
		if (x_ce_max < x_pos_end - 1) ClusterSubSeq(patient_info, seq2d, x_ce_max, x_pos_end, y_cs_max, y_ce_max);
		if (x_cs_max > x_pos_start +1 && y_ce_max < y_pos_end - 1) ClusterSubSeq(patient_info, seq2d, x_pos_start, x_cs_max, y_ce_max, y_pos_end);
		if (y_ce_max < y_pos_end - 1) ClusterSubSeq(patient_info, seq2d, x_cs_max, x_ce_max, y_cs_max, y_ce_max);
		if (x_ce_max < x_pos_end - 1 && y_ce_max < y_pos_end - 1) ClusterSubSeq(patient_info, seq2d, x_ce_max, x_pos_end, y_ce_max, y_pos_end);


 	}

	//cout << "bp6" << endl;

	return 1;
}




//Model averaging based on all possible candidate models
int Cluster::ModelAveraging(long x_pos_start, long x_pos_end, long y_pos_start, long y_pos_end, long x_cs_max, long x_ce_max, long y_cs_max, long y_ce_max, double p0_max, double pc_max, double MIN_cri) {
	//cout << "b1" << endl;
	long i, j, k, l;
	char spot = 'H';
	if(pc_max<p0_max) spot = 'C';



	//for (i=pos_start; i<=pos_end; i++) vec_MS_rate[i] = p0_max;



	for (i = x_pos_start; i <= x_pos_end; i++) {
		for (j = y_pos_start; j <= y_pos_end; j++) {
			vec_MS_rate[j][i] = p0_max;
		}
	}

	for (i=x_cs_max; i<=x_ce_max; i++) {
		for (j = y_cs_max; j <= y_ce_max; j++) {
			vec_MS_rate[j][i] = pc_max;
			vec_spot[j][i] = spot;
		}
	}

	//cout << "b2" << endl;

	/* Model Averaging */
	double all_weight = 0.0;
	for (i=0; i<vec_AllModels.size(); i++) {		
		vec_AllModels[i].CW = exp(-0.5*(vec_AllModels[i].CW - MIN_cri));
		all_weight += vec_AllModels[i].CW;
	}

	for (i=0; i<vec_AllModels.size(); i++) {		
		vec_AllModels[i].CW = vec_AllModels[i].CW / all_weight;
	}
	


	cout << "CI length is " << endl;


	vector <double> weights;



	cout << vec_AllModels.size() << endl;



	for (j=0; j<vec_AllModels.size(); j++) {
		//cout << vec_AllModels[j].CW << endl;
		weights.push_back(vec_AllModels[j].CW);
	}		

	//stable_sort(weights.begin(), weights.end());


	double sum_of_weights = 0.0;

	for (j = 0; j < weights.size(); ++j) {
		//cout << weights[j] << endl;
		sum_of_weights += weights[j];
	}



	vector <int> weight_indices_sorted;

	//std::vector<int> y(x.size());
	
	for (j = 0; j < weights.size(); ++j) {
		weight_indices_sorted.push_back((int) j);
	}

	//auto comparator = (int a, int b){ return weights[a] < weights[b]; };
	sort(weight_indices_sorted.begin(), weight_indices_sorted.end(), [&weights](int a, int b) { return weights[a] < weights[b]; });


	//cout << sum_of_weights << endl;


	//for (j = 0; j < weights.size(); ++j) {
		//cout << weight_indices_sorted[j] << endl;
		//sum_of_weights += weights[j];
	//}	

	vector <double> sorted_weights;
	for (j = 0; j < weights.size(); ++j) {
		sorted_weights.push_back(weights[weight_indices_sorted[j]]);
		//cout << sorted_weights[j] << endl;
	}


	vector <CandidateModels> sampled_models;

	srand((unsigned) time(NULL));

	double cumulative;

	double random_number;
	int index;
	for (j = 0; j < sample_size; ++j) {

		random_number = ((double) rand())/(double) RAND_MAX ;
		//cout << random_number << endl;

		cumulative = 0.0;
		index = -1;

		while (cumulative < random_number) {
			++index;
			cumulative += sorted_weights[weights.size() - 1 - index];


		}
		//cout << index << endl;
		sampled_models.push_back(vec_AllModels[weight_indices_sorted[weights.size() -1 - index]]);


	}

	double sampled_models_weight_sum = 0.0;

	for (j = 0; j < sampled_models.size(); ++j) {
		sampled_models_weight_sum += sampled_models[j].CW;


	}


	for (j = 0; j < sampled_models.size(); ++j) {
		sampled_models[j].CW  = sampled_models[j].CW/sampled_models_weight_sum;

		
	}
	



	for (i=x_pos_start; i<=x_pos_end; i++) {
		for (k = y_pos_start; k <= y_pos_end; k++) {
			double rate = 0.0;
			vector<CI> CIs;//probability and weight for all possible models
			vector <CI> CIs2;
			//cout << "b4" << endl;

			for (j=0; j<sampled_models.size(); j++) {

				double site_weight = sampled_models[j].CW;
				double site_rate = sampled_models[j].pc;
				if (i<sampled_models[j].x_cs || i>sampled_models[j].x_ce || k < sampled_models[j].y_cs || k > sampled_models[j].y_ce) {
					site_rate = sampled_models[j].p0;
				}
				
				CI tmp_cis(site_weight, site_rate);
				CIs.push_back(tmp_cis);
				CIs2.push_back(tmp_cis);

				rate += site_rate*site_weight;
			}		
			vec_MA_rate[k][i] = rate;




			vec_lower_rate[k][i] = 0.0;
			vec_upper_rate[k][i] = 1.0;
			//sort by weight
			double lower=0.0, upper=0.0;
			  int flag_lower=0;
			  //Here confidece_interval=0.025
			 //cout << "b41" << endl;
				while (lower<confidence_interval) {
					//cout << "b42" << endl;
					long min_pos;//, max_pos;
					min_pos = -1;//max_pos = -1;
					j = 0;
					long min = 2;
					//long max = -1;
					while (j<CIs.size()) {
						if(CIs[j].p<min && CIs[j].p != -1) {
							min_pos = j;
							min = CIs[j].p;
						}
						// if(CIs[j].p>max && CIs[j].p != -1) {
						// 	max_pos = j;
						// 	max = CIs[j].p;
						// }
						j++;
					}
					if (min_pos == -1) {// && max_pos == -1) {
						break;
					}


					//cout << "b43" << endl;
					//if (lower <confidence_interval) {
					  // if(max_pos>min_pos){
					  //   flag_lower=1;
					  // }
						lower += CIs[min_pos].weight;
						//cout << "b431" << endl;
						//fflush(stdout);
						//vec_lower_rate[k][i] = CIs[min_pos].p;
						if (CIs[min_pos].p <= vec_MA_rate[k][i]) {
							vec_lower_rate[k][i] = CIs[min_pos].p;
						}
						//cout << "b432" << endl;
						//cout << CIs.size() << endl;
						//fflush(stdout);
						//CIs.erase(CIs.begin() + min_pos);
						CIs[min_pos].p = -1;
						//cout << "b433" << endl;
						//fflush(stdout);
					//}
					//cout << "b44" << endl;
					// if(upper <confidence_interval) {
					//   // if(flag_lower==1 && max_pos!=0){
					//   //   max_pos--;
					//   // }
					// 	upper += CIs[max_pos].weight;
					// 	vec_upper_rate[k][i] = CIs[max_pos].p;
					// 	//CIs.erase(CIs.begin() + max_pos);
					// 	CIs[max_pos].p = -1;
					// }
					//flag_lower=0;
				}

				while (upper<confidence_interval) {
					//cout << "b42" << endl;
					long max_pos;//, max_pos;
					max_pos = -1;//max_pos = -1;
					j = 0;
					long max = -1;
					//long max = -1;
					while (j<CIs2.size()) {
						if(CIs2[j].p>max && CIs2[j].p != -1) {
							max_pos = j;
							max= CIs2[j].p;
						}
						// if(CIs[j].p>max && CIs[j].p != -1) {
						// 	max_pos = j;
						// 	max = CIs[j].p;
						// }
						j++;
					}
					if (max_pos == -1) {// && max_pos == -1) {
						break;
					}


					//cout << "b43" << endl;
					//if (lower <confidence_interval) {
					  // if(max_pos>min_pos){
					  //   flag_lower=1;
					  // }
						upper += CIs2[max_pos].weight;
						//cout << "b431" << endl;
						//fflush(stdout);
						if (CIs2[max_pos].p >= vec_MA_rate[k][i]) {
							vec_upper_rate[k][i] = CIs2[max_pos].p;
						}
						//cout << "b432" << endl;
						//cout << CIs.size() << endl;
						//fflush(stdout);
						//CIs.erase(CIs.begin() + min_pos);
						CIs2[max_pos].p = -1;
						//cout << "b433" << endl;
						//fflush(stdout);
			
				}
			//   cout << "ABOUT TO PRINT" << endl;
			// cout << CIs.size() << endl;
			// for (l = 0; l < CIs.size(); ++l) {
			// 	cout << CIs[l].p << "\t" << CIs[l].weight;
			// }



			//cout << "b45" << endl;		
		

		CIs.clear();
		}
	}

	cout << "b5" << endl;

	return 1;
}


/**************************************************
* Function: parseParameter
* Input Parameter: int, const char* []
* Output: Parse the input parameters
* Return Value: bool 
***************************************************/
bool Cluster::parseParameter(int argc, const char* argv[]) {
	
	bool flag = true;
	
	try {
		//Input sequence
		int inputfile_flag=0;
		
		if (argc!=3 && argc!=5 && argc!=7 && argc!=9) {
			throw 1;			
		}		
		else 
		{
			int i;
			string temp;			
			//parse parameters			
			for (i=1; i<argc; i++) {				
				temp = stringtoUpper(argv[i]);
				//Input fasta file
				if (temp=="-I" && (i+1)<argc && inputfile_flag==0) {
					seqfile = argv[++i];						
					inputfile_flag++;
				}

				//BIC/AIC/AICc
				else if (temp=="-C" && (i+1)<argc) {
					int num = CONVERT<int>(argv[++i]);
					if (num==0 || num==1 || num==2) {
						criterion_type = num;
					}
					else {
						throw 1;
					}
				}				
				//Sequence type
				else if (temp=="-T" && (i+1)<argc) {
					int num = CONVERT<int>(argv[++i]);
					if (num==0 || num==1) {
						seq_type = num;
					}
					else {
						throw 1;
					}
				}
				else if (temp == "-S" && (i+1) < argc) {
					int num = CONVERT<int>(argv[++i]);
					sample_size = num;
				}
				//Model selection and model averaging
				else if (temp=="-M" && (i+1)<argc) {
					int num = CONVERT<int>(argv[++i]);
					if (num==0 || num==1) {
						MS_only = num;
					}
					else {
						throw 1;
					}
				}
				else {
					throw 1;
				}
			}			
		}
		//no input file
		if (inputfile_flag==0) throw 1;
	}
	catch (...) {
		flag = false;		
	}
	
	return flag;
}
