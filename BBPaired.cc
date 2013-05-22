#include "pedParser.h"
#include "BBPaired.h"
#include <getopt.h>
#include <fstream>


BBPaired::BBPaired(int argc, char* argv[])
{ 
  ip_RD_f = "empty";
  op_calls_f = "empty";
  ped_f = "empty";
  alphaRR = -1;
  betaRR = -1;
  alphaRA = -1;
  betaRA = -1;
  alphaAA = -1;
  betaAA = -1;
  parseArguments(argc, argv);
  getSampleIDFromPED();
  parseFile();
}


int BBPaired::parseArguments(int argc, char* argv[])
{
  static struct option long_options[] = {
                                          {"input_RD_file", 1, 0, 0}, 
					  {"op_calls_file", 1, 0, 1},
					  {"ped", 1, 0, 2},
					  {"alphaRR", 1, 0, 3},
					  {"betaRR", 1, 0, 4},
					  {"alphaRA", 1, 0, 5},
					  {"betaRA", 1, 0, 6},
					  {"alphaAA", 1, 0, 7},
					  {"betaAA", 1, 0, 8}
  };
  optind = 0;
  // Read in Command Line arguments
  while (1) {
    int option_index = 0;
    int c = getopt_long (argc, argv, "", long_options, &option_index);
    if (c == -1)
      break;
    switch(c) 
      {      
        case 0:
	  ip_RD_f = optarg;
	  break;

        case 1:
	  op_calls_f = optarg;
	  break;

        case 2:
	  ped_f = optarg;
	  break;

        case 3:
	  alphaRR = atof(optarg);
	  break;

        case 4:
	  betaRR = atof(optarg);
	  break;

        case 5:
	  alphaRA = atof(optarg);
	  break;

        case 6:
	  betaRA = atof(optarg);
	  break;

        case 7:
	  alphaAA = atof(optarg);
	  break;

        case 8:
	  betaAA = atof(optarg);
	  break;

        default:
	  break;
      }
  }
  return 0;
}

int BBPaired::getSampleIDFromPED()
{
  Trio* trios;
  Pair* pairs;
  int trio_count = 0, pair_count = 0;
  parse_ped(ped_f, &trios, &pairs, trio_count, pair_count);
  tumorID = pairs[0].tumorID;
  normalID = pairs[0].normalID;
}


int BBPaired::parseFile()
{
  ifstream fin(ip_RD_f.c_str());
  if(!fin) {
    cerr<<"Unable to open file: "<<ip_RD_f<<".Exiting !"<<endl;
    exit (EXIT_FAILURE);
  }
  string line, sampleID;
  while(fin.good()) {
    getline(fin, line);
    sampleID = parseLine(line);
    if(sampleID == tumorID) // do it for just the first pair for now.
      calcBBLik(tumor_pt);
    else if(sampleID == normalID) // do it for just the first pair for now.
      calcBBLik(normal_pt);
  }
  return 0;
}


string BBPaired::parseLine(string line) //line is of the format "chr\tpos\tref\tref_RD\talt\talt_RD"
{
  string sampleID;
  VI1.reset();
  istringstream iss(line);
  iss >> VI1.chr;
  iss >> VI1.pos;
  iss >> VI1.ref;
  iss >> VI1.ref_RD;
  iss >> VI1.alt;
  iss >> VI1.alt_RD;  
  iss >> sampleID;
  return sampleID;
}


int BBPaired::calcBBLik(pair_t& pt)
{
  VI1.print();
  return 0;
}

