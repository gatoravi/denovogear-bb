#ifndef BB_PAIRED_H
#define BB_PAIRED_H

#include<iostream>
#include "denovogear.h"
#include "makeLookup.h"
#include "pedParser.h"
using namespace std;

struct variant_info{
  string chr;
  long pos;
  string ref;
  string alt;
  int ref_RD; // no of reads supporting ref
  int alt_RD; // no of reads supporting alt 
  void reset()
  {
    chr = "NA";
    pos = -1;
    ref = "NA";
    alt = "NA";
    ref_RD = -1;
    alt_RD = -1;
  }
  void print()
  {
    cout<<"\nV1\t"<<chr<<"\t"<<pos<<"\t"<<ref<<"\t"<<ref_RD;
    cout<<"\t"<<alt<<"\t"<<alt_RD<<endl;
  }

};


class BBPaired {
  public:
    BBPaired(int argc, char* argv[]);
  private:
    string ip_RD_f;
    string op_calls_f;
    string ped_f;
    string tumorID;
    string normalID;
    double alphaRR;
    double betaRR;
    double alphaRA;
    double betaRA;
    double alphaAA;
    double betaAA;
    pair_t tumor_pt;
    pair_t normal_pt;
    variant_info VI1;
    Pair* pairs;
    int pair_count;
    int getSampleIDFromPED();
    int parseArguments(int argc, char* argv[]);
    int calcBBLik(pair_t& pt);
    int parseFile();
    string parseLine(string line);
};
#endif
