#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <armadillo>
#include <unistd.h>
#include <vector>
#include "util.h"
#include "postcal.h"
#include "model.h"
#include "PIPSORTConfig.h"

using namespace std;

/*
 Reads the content of the file, return a vector of paths
 @param fileName the name of file that contains all the paths to z or ld files
 @return vector of paths
 */
vector<string> read_dir(string fileName){
    vector<string> dirs;
    string data;

    ifstream fin(fileName.c_str(), std::ifstream::in);
    if (!fin) {
        cerr << "Error: unable to open " << fileName << endl;
        exit(1); // terminate with error
    }

    while(fin.good()){
        getline(fin,data);
        if(data != "") {
            dirs.push_back(data); }}
    fin.close();
    return dirs;
}

vector<int> read_sigma(string sample_size) {
    vector<int> sizes;
    string current_size = "";
    for (int i=0; i < sample_size.size(); i++) {
        if (sample_size[i] == ',') {
            sizes.push_back(stod(current_size)); // push back previous sample size
            current_size = ""; // reset
        }
        else if (isdigit(sample_size[i])) {
            current_size += sample_size[i];
        }
        else {
            cerr << "Error: sample size is not in the right format" << endl;
            exit(1);
        }
    }
    if (current_size != "") {
        sizes.push_back(stod(current_size)); // push back last sample size
    }
    return sizes;
}

void show_help() {
    std::stringstream help_msg;
    help_msg << "\nUsage: PIPSORT [OPTIONS] "
             << "-l <LDFILE> "
             << "-z <ZFile> "
             << "-m <snpMapFile> "
             << "-n <int,int> "
             << "-o <string>"
             << "\n\n Required options:\n"
             << "\t" << "-l         <LDFILE>     " << "\t" << "File containing paths to ld files" << "\n"
             << "\t" << "-z         <ZFILE>      " << "\t" << "File containing paths to Z files" << "\n"
             << "\t" << "-m         <snpMapFile> " << "\t" << "File mapping indexes to SNPs in each study" << "\n"
             << "\t" << "-n         <int1,int2>  " << "\t" << "Sample sizes (integers) of each study. e.g. 50,100" << "\n"
             << "\t" << "-o         <OUTPREFIX>  " << "\t" << "Prefix for output files" << "\n"
             << "\n Additional optional parameters:\n"
             << "\t" << "-p         <SHAREPARAM> " << "\t" << "Sharing parameter (default 0.75)" << "\n"
             << "\t" << "-c         <NUMCAUSAL>  " << "\t" << "set the maximum number of causal SNPs (default 3)" << "\n"
             //<< "\t" << "-k         <NUMCAUSAL>  " << "\t" << "set the number of causal SNPs per study (default 3)" << "\n"
             << "\t" << "-r         <RHO>        " << "\t" << "set $rho$ probability (default 0.95)" << "\n"
             << "\t" << "-g         <GAMMA>      " << "\t" << "set $gamma$ the prior of a SNP being causal (default 0.01)" << "\n"
             << "\t" << "-t         <TAU_SQR>    " << "\t" << "set the heterogeneity (t^2) across studies, (default 0.52)" << "\n"
             << "\t" << "-s         <SIGMA_GSQR> " << "\t" << "set the NCP variance for the smallest study, (default 5.2)" << "\n"
             << "\t" << "-a         <THRESHOLD>  " << "\t" << "if a variant has a posterior below this threshold, " << "\n"
             << "\t" << "                        " << "\t" << "do not include it in the causal set (default 0)" << "\n"
             //<< "\t" << "-f         <int>        " << "\t" << "to out the probaility of different number of causal SNP" << "\n"
             << "\n Options for stochastic shotgun search (SSS):\n"
             << "\t" << "-q         <INT>        " << "\t" << "set to 1 to perform SSS" << "\n"
             << "\n Options for adding custom configurations to test:\n"
             << "\t" << "-b         <CONFIGFILE> " << "\t" << "Optional causal configuration file" << "\n"
             << "\t" << "-d         <NUMCONFIG>  " << "\t" << "Number of configurations (rows in configuration file)" << "\n"
             << "\t" << "-e         <NUMGROUPS>  " << "\t" << "Number of groups (columns in configuration file)" << "\n"
             << "\n Additional options:\n"
             << "\t" << "-h                      " << "\t" << "display this help screen" << "\n"
             << "\t" << "-v                      " << "\t" << "show the version number" << "\n";

    cerr << help_msg.str();
    exit(0);
}

int main( int argc, char *argv[]  ){
    int totalCausalSNP = 3;
    double gamma = 0.01;
    double sharing_param = 0.75;
    double rho = 0.95;
    bool histFlag = false;
    int oc = 0;
    double tau_sqr = 0.52;
    double sigma_g_squared = 5.2;
    double cutoff_threshold = 0;

    string ldFile = "";
    string zFile  = "";
    string snpMapFile = "";
    string outputFileName = "";
    string sample_s = "";
    string num_causal_s = "";
    string configsFile = "";
    int num_groups = 0; //num columns in configsFile
    int num_configs = 0; //num rows in configsFile
    int sss_flag = 0;

    while ((oc = getopt(argc, argv, "vhl:o:z:m:p:r:c:k:g:f:t:s:n:a:b:d:e:q:x")) != -1) {
	    //TODO P3 last char in this colon separated list does not work, optarg comes in as null. -x is dummy flag. Should it be :x: (colon at end)?
        if (oc != 'v' && oc != 'h' && oc != 'x') {
	       if ( (optarg == NULL) || (*optarg == '\0') ) {
               printf("optarg is NULL\n");
	           exit(1);
	       }
        }
        switch (oc) {
            case 'v':
                cerr << PIPSORT_VER << endl;
		        exit(0);
            case 'h':
                show_help();
                exit(0);
            case 'l':
                ldFile = string(optarg);
                break;
            case 'o':
                outputFileName = string(optarg);
                break;
            case 'z':
                zFile = string(optarg);
                break;
            case 'm':
		        snpMapFile = string(optarg);
            case 'n':
                sample_s = string(optarg);
                break;
	        // optional argument: file with configuration vectors
	        case 'b':
	            configsFile = string(optarg);
	            break;
	        case 'd':
	            num_configs = atoi(optarg);
	            break;
	        case 'e':
	            num_groups = atoi(optarg);
	            break;
            // optional arguments: parameters for fine mapping
	        case 'p':
		        sharing_param = atof(optarg);
		        break;
            case 'r':
                rho = atof(optarg);
                break;
            case 'c':
                totalCausalSNP = atoi(optarg);
                break;
            case 'k':
                num_causal_s = string(optarg);
                break;
            case 'g':
                gamma = atof(optarg);
                break;
            case 'f':
                histFlag = true;
                break;
            case 't':
                tau_sqr = atof(optarg);
                break;
            case 's':
                sigma_g_squared = atof(optarg);
                break;
	        case 'q':
		        sss_flag = stoi(optarg);
		        break;
            case ':':
            case '?':
            case 'a':
                cutoff_threshold = atof(optarg);
                break;
            default:
                break;
        }
    }

    if (ldFile == "" or zFile == "" or snpMapFile == "" or outputFileName == "" or sample_s == "") {
        cerr << "Error: -l, -z, -o, and -n are required" << endl;
        show_help();
        exit(1);
    }

    if ( configsFile != "" ) {
      if (num_configs <= 0) {
        cerr << "Number of configs must be greater than 0" << endl;
	    exit(1);
      }
      if (num_groups <= 0) {
        cerr << "Number of groups must be greater than 0" << endl;
        exit(1);
      }
    }

    bool do_sss = false;
    if ( sss_flag == 1 ) {
       do_sss = true;
    }

    vector<string> ldDir = read_dir(ldFile);
    vector<string> zDir = read_dir(zFile);
    vector<int> sample_sizes = read_sigma(sample_s);
    vector<int> num_causal;
    static const int finalTotalCausalSNP = totalCausalSNP;

    if (num_causal_s != "") {
	for (int i = 0; i < sample_sizes.size(); i++) {
	    num_causal.push_back(3);
	}
    } else {
        num_causal = read_sigma(num_causal_s);
    }

    if (ldDir.size() != zDir.size() || ldDir.size() != sample_sizes.size()) {
        cerr << "Error: LD files, Z files, and sample sizes do not match in number" << endl;
        exit(1);
    }  
    omp_set_num_threads(1);

    Model mpipsort(ldDir, zDir, snpMapFile, configsFile, num_configs, num_groups, do_sss, sample_sizes, num_causal, outputFileName, finalTotalCausalSNP, sharing_param, rho, histFlag, gamma, tau_sqr, sigma_g_squared, cutoff_threshold);
    mpipsort.run();
    mpipsort.finishUp();
    return 0;
}
