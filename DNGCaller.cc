#include "DNGCaller.h"
#include "common.h"
#include "lookup.h"
#include "denovogear.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include <math.h>
#include <cstring>
#define MIN_PHRED_LIKE 255
#define MAX_PHRED_LIKE 0

DNGCaller::DNGCaller(string tumor_mp, string normal_mp, 
		     std::vector<double> tumor_ab1, 
		     std::vector<double> normal_ab1)
{
  //cout<<"constructing ";
  //cout<<"\n1t size "<<tumor_ab1.size()<<"\n";
  tumor_mp_f = tumor_mp;
  normal_mp_f = normal_mp;
  tumor_ab = tumor_ab1;
  normal_ab = normal_ab1;
  n_site_pass = 0;
  callMakePairedLookup(tgtPair1, lookupPair1);
  iterateMpileups();
}

DNGCaller::~DNGCaller()
{
}

int DNGCaller::iterateMpileups()
{
  //cout<<"\nt size "<<tumor_ab.size()<<"\n";
  ifstream tin1(tumor_mp_f.c_str());
  if(!tin1) {
    cerr<<"Unable to open "<<tumor_mp_f<<" .Exiting !";
    exit(EXIT_FAILURE);
  }
  ifstream nin1(normal_mp_f.c_str());
  if(!nin1) {
    cerr<<"Unable to open "<<normal_mp_f<<" .Exiting !";
    exit(EXIT_FAILURE);
  }
  string l_t, l_n;
  string t_chr, n_chr;
  long t_pos, n_pos;
  getline(tin1, l_t);
  getline(nin1, l_n);

  while(tin1.good() && nin1.good()) {
    istringstream iss_t(l_t);
    istringstream iss_n(l_n);
    iss_t >> t_chr;
    iss_t >> t_pos;
    iss_n >> n_chr;
    iss_n >> n_pos;
    if(t_chr == n_chr && t_pos == n_pos) {
      calcLik(l_t, l_n);
      getline(tin1, l_t);
      getline(nin1, l_n);
    }
    else if(t_chr == n_chr) {
	if(t_pos > n_pos) {
	  getline(nin1, l_n);
	  continue;
	}
	else {
	  getline(tin1, l_t);
	  continue;
	} 
    }
    else if(t_chr > n_chr) {
      getline(nin1, l_n);
      continue;
    }
    else if(n_chr > t_chr){
      getline(tin1, l_t);
      continue;
    } 
  }
  return 0;
}

int DNGCaller::parseLine(string l, t_variant_info& vi1)
{
  cout<<l<<endl;
  istringstream iss(l);
  string chr, ref, pileup, field;
  long pos;
  int RD, a_count, c_count, g_count, t_count, ref_count;
  typedef std::pair<std::string, int> basecount;
  std::vector<basecount> basecounts;

  iss >> chr;
  iss >> pos;
  iss >> field;
  RD = atoi(field.c_str());
  iss >> ref;
  iss >> pileup;
  iss >> field;
  a_count = atoi(field.c_str());
  iss >> field;
  c_count = atoi(field.c_str());
  iss >> field;
  g_count = atoi(field.c_str());
  iss >> field;
  t_count = atoi(field.c_str());

  if(ref == "a" || ref == "A") {
    ref = "A";
    ref_count = a_count;
  }
  if(ref == "c") {
    ref = "C";
    ref_count = c_count;
  }
  if(ref == "g") {
    ref = "G";
    ref_count = g_count;
  }
  if(ref == "t") {
    ref = "T";
    ref_count = t_count;
  }

  basecounts.push_back(std::make_pair("A", a_count));
  basecounts.push_back(std::make_pair("C", c_count));
  basecounts.push_back(std::make_pair("G", g_count));
  basecounts.push_back(std::make_pair("T", t_count));

  std::sort (basecounts.begin(), basecounts.end(), compareBaseCounts);
  /*std::cout<<endl<<"Line "<<l;
  std::cout<<endl<<basecounts[0].first<<"\t"<<basecounts[0].second;
  std::cout<<endl<<basecounts[1].first<<"\t"<<basecounts[1].second;
  std::cout<<endl<<basecounts[2].first<<"\t"<<basecounts[2].second;
  std::cout<<endl<<basecounts[3].first<<"\t"<<basecounts[3].second;*/
  
  vi1.chr = chr;
  vi1.pos = pos;
  vi1.ref = ref;
  vi1.ref_RD = ref_count;
  if(basecounts[0].first == ref) {
    vi1.alt = basecounts[1].first;
    vi1.alt_RD = basecounts[1].second;
  }
  else {
    vi1.alt = basecounts[0].first;
    vi1.alt_RD = basecounts[0].second;
  }
  //vi1.print();
  return 0;
}

double DNGCaller::lbetaFn(double alpha, double beta)
{
  return lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
}

double DNGCaller::lik_bb(double alpha, double beta, int n, int k)
{
  while(k>500 || n>500) {
    k = k/2;
    n = n/2;
  }
  double t1  = k + alpha;
  double t2  = n - k + beta;
  double lnum = lbetaFn(t1, t2);		
  double ldenom = lbetaFn(alpha, beta);
  double log_lik  = lnum - ldenom;
  return log_lik;
}

int DNGCaller::convertVI2paired(t_variant_info& vi1, pair_t& pt1, std::vector<double>& ab, string sampleID)
{
  std::map<string, int> GT_lik;
  GT_lik["AA"] =  MIN_PHRED_LIKE;
  GT_lik["AC"] =  MIN_PHRED_LIKE;
  GT_lik["AG"] =  MIN_PHRED_LIKE;
  GT_lik["AT"] =  MIN_PHRED_LIKE;
  GT_lik["CC"] =  MIN_PHRED_LIKE;
  GT_lik["CG"] =  MIN_PHRED_LIKE;
  GT_lik["CT"] =  MIN_PHRED_LIKE;
  GT_lik["GG"] =  MIN_PHRED_LIKE;
  GT_lik["GT"] =  MIN_PHRED_LIKE;
  GT_lik["TT"] =  MIN_PHRED_LIKE;

  string RR = vi1.ref + vi1.ref;
  string RA = vi1.ref + vi1.alt;
  string AA = vi1.alt + vi1.alt;

  //std::cout<<"RR "<<RR;
  //std::cout<<"RA "<<RA;
  //std::cout<<"AA "<<AA;

  //std::cout<<"1"<<endl;

  //std::cout<<ab.size()<<endl;

  double log_lik_bb_RR   =  -10 * lik_bb(ab[0], ab[1], vi1.ref_RD + vi1.alt_RD, vi1.alt_RD);
  double log_lik_bb_RA   =  -10 * lik_bb(ab[2], ab[3], vi1.ref_RD + vi1.alt_RD, vi1.alt_RD); 
  double log_lik_bb_AA   =  -10	* lik_bb(ab[4], ab[5], vi1.ref_RD + vi1.alt_RD, vi1.alt_RD); 

  cout<<"\nbefore norm liks "<<log_lik_bb_RR<<"\t"<<log_lik_bb_RA<<"\t"<<log_lik_bb_AA;  
  double norm_bb_liks[3];
  normalize_liks(log_lik_bb_RR, log_lik_bb_RA, log_lik_bb_AA, norm_bb_liks);
  cout<<"\nafter norm liks "<<norm_bb_liks[0]<<"\t"<<norm_bb_liks[1]<<"\t"<<norm_bb_liks[2];

  GT_lik[RR] = norm_bb_liks[0];
  GT_lik[RA] = norm_bb_liks[1];
  GT_lik[AA] = norm_bb_liks[2];
  
  strcpy(pt1.chr, vi1.chr.c_str());
  pt1.pos = vi1.pos;
  pt1.ref_base = vi1.ref[0];
  strcpy(pt1.alt, vi1.alt.c_str());
  pt1.depth = vi1.ref_RD + vi1.alt_RD;
  pt1.rms_mapQ = -1;
  strcpy(pt1.id, sampleID.c_str());
  std::map<std::string, int>::iterator it;
  int index = 0;
  for(it = GT_lik.begin(); it != GT_lik.end(); it++) {
    cout<<"\nIT "<<it->first<<"\t"<<it->second;
    pt1.lk[index++] = it->second;
  }
  cout<<endl;
  
  return 0;
}

int DNGCaller::normalize_liks(double log_lik_RR, double log_lik_RA, double log_lik_AA, double* normalized_liks)
{
  double min_log_lik = min(log_lik_RR, min(log_lik_RA, log_lik_AA));
  double norm_log_lik_RR = (int)(log_lik_RR - min_log_lik + 0.499);
  double norm_log_lik_RA = (int)(log_lik_RA - min_log_lik + 0.499);
  double norm_log_lik_AA = (int)(log_lik_AA - min_log_lik + 0.499);
  			
  if(norm_log_lik_RR > 255)
    norm_log_lik_RR = MIN_PHRED_LIKE;

  if(norm_log_lik_RA > 255)
    norm_log_lik_RA = MIN_PHRED_LIKE;

  if(norm_log_lik_AA > 255)
    norm_log_lik_AA = MIN_PHRED_LIKE;

  if(norm_log_lik_RR < 0)
    norm_log_lik_RR = MAX_PHRED_LIKE;

  if(norm_log_lik_RA < 0)
    norm_log_lik_RA = MAX_PHRED_LIKE;

  if(norm_log_lik_AA < 0)
    norm_log_lik_AA = MAX_PHRED_LIKE;

  normalized_liks[0] = norm_log_lik_RR;
  normalized_liks[1] = norm_log_lik_RA;
  normalized_liks[2] = norm_log_lik_AA;

  return 0;
}


int DNGCaller::calcLik(string l_t, string l_n)
{
  //cout<<"\ntumor  "<<l_t;
  //cout<<"\nnormal "<<l_n;
  t_variant_info tvi, nvi;
  pair_t pt_t, pt_n;
  cout<<endl<<"tumor";
  parseLine(l_t, tvi);
  convertVI2paired(tvi, pt_t, tumor_ab, "tumor");
  cout<<endl<<"normal";
  parseLine(l_n, nvi);
  convertVI2paired(nvi, pt_n, normal_ab, "normal");
  double pp_cutoff1 = 0.001;
  int RD_cutoff1= 10;
  int flag1 = 0;
  ofstream fo_vcf1;
  string op_vcf_f1 = "empty";
  // Create paired lookup
  pair_like(pt_t, pt_n, tgtPair1, lookupPair1, flag1, op_vcf_f1, fo_vcf1, pp_cutoff1, RD_cutoff1, n_site_pass);
  return 0;
}
