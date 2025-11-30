#ifndef POSTCAL_H
#define POSTCAL_H

#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <unordered_map>

#include <armadillo>

using namespace std;
using namespace arma;

void printGSLPrint(mat A, int row, int col);

class PostCal{
private:

    double gamma;        // the probability of SNP being causal
    int totalSnpCount;
    double * postValues;    //the posterior value for each SNP being causal
    double * noCausal;
    double * sharedPips;
    double * sharedLL;
    double * notSharedLL;
    double * histValues;    //the probability of the number of causal SNPs, we make the histogram of the causal SNPs
    int snpCount;        //total number of variants (SNP) in a locus
    string configsFile; //optional configs file, will be empty string if not provided
    int num_configs;
    int num_groups;
    const int maxCausalSNP;    //maximum number of causal variants to consider in a locus
    double sigmaDet;    //determinant of matrix

    double totalLikeLihoodLOG; //Compute the total log likelihood of all causal status (by likelihood we use prior)
    double sharing_param;
    double t_squared;    //tau^2 (heterogeneity)
    double s_squared;    //sigma_g^2 (heritability)
    const int num_of_studies;
    bool haslowrank = false;

    mat sigmaMatrix;
    mat invSigmaMatrix;
    mat statMatrix;
    mat statMatrixtTran;
    vector<vector<string> > * SNP_NAME;
    vector<int> sample_sizes;
    vector<int> num_snps_all;
    vector<int> num_causal;
    vector<vector<int>> idx_to_snp_map;
    vector<vector<int>> idx_to_union_pos_map;
    vector<string> all_snp_pos;
    int unionSnpCount;
    std::unordered_map<vector<int>, double, VecHash> config_hashmap;

    //addition in log space
    double addlogSpace(double a, double b) {
        if (a == 0)
            return b;
        if (b == 0)
            return a;
        double base = max(a,b);
        if (base - min(a,b) > 700)
            return base;
        
        return(base + log(1+exp(min(a,b)-base)));
    }

public:
    /*
     constructor
    */
    PostCal(mat * BIG_SIGMA, vector<double> * S_LONG_VEC, int snpCount, string configsFile, int num_configs, int num_groups, const int MAX_causal, vector<int> num_causal, vector<vector<string> > * SNP_NAME, double sharing_param, double gamma, double t_squared, double s_squared, const int num_of_studies, vector<int> sample_sizes, vector<int> num_snps_all, bool lowrank, vector<vector<int>> idx_to_snp_map, vector<vector<int>> idx_to_union_pos_map, vector<string> all_snp_pos) : maxCausalSNP(MAX_causal),num_of_studies(num_of_studies){
        this->gamma = gamma;
        this->SNP_NAME = SNP_NAME;
        this-> snpCount = snpCount;
	this->configsFile = configsFile;
	this->num_configs = num_configs;
	this->num_groups = num_groups;
	this-> totalSnpCount = std::accumulate(num_snps_all.begin(), num_snps_all.end(), 0);
        //this-> maxCausalSNP = MAX_causal;
        //this-> postValues = new double [snpCount];
        this-> postValues = new double [totalSnpCount];
        for(int i = 0; i < totalSnpCount; i++)
            this->postValues[i] = 0;
	this-> noCausal = new double [num_of_studies];
	for ( int i = 0; i < num_of_studies; i++ ) {
           this->noCausal[i] = 0;
	}
        this-> histValues = new double [MAX_causal+1];
        for(int i= 0; i <= maxCausalSNP;i++)
            this->histValues[i] = 0;
        this-> totalLikeLihoodLOG = 0;
	this-> sharing_param = sharing_param;
        this-> t_squared = t_squared;
        this-> s_squared = s_squared;
        //this-> num_of_studies = num_of_studies;
        this-> sample_sizes = sample_sizes;
	this-> num_causal = num_causal;
	this-> num_snps_all = num_snps_all;
        this-> haslowrank = lowrank;
	this-> idx_to_snp_map = idx_to_snp_map;
	this-> idx_to_union_pos_map = idx_to_union_pos_map;
	this-> all_snp_pos = all_snp_pos;
	this-> unionSnpCount = all_snp_pos.size();
	this-> sharedPips = new double [unionSnpCount];
	this-> sharedLL = new double [unionSnpCount];
	this-> notSharedLL = new double [unionSnpCount];
	for ( int i = 0; i < unionSnpCount; i++ ) {
          this->sharedPips[i] = 0;
	  this->sharedLL[i] = 0;
	  this->notSharedLL[i] = 0;
	}

        // statMatrix is the z-score matrix of mn*1, m = number of snps, n = num of studies
        statMatrix = mat (totalSnpCount, 1);
        statMatrixtTran = mat (1, totalSnpCount);
        for(int i = 0; i < totalSnpCount; i++) {
            statMatrix(i,0) = (*S_LONG_VEC)[i];
            statMatrixtTran(0,i) = (*S_LONG_VEC)[i];
        }
        // sigmaMatrix is diagonal matrix of sigma matrices for each study i, same for invSigmaMatrix, sigmaDet
        sigmaMatrix = mat (totalSnpCount, totalSnpCount);
        sigmaMatrix = (*BIG_SIGMA);
        /*
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0, 1);
        for(int i = 0; i < snpCount * num_of_studies; i++) {
            for (int j = 0; j < snpCount * num_of_studies; j++) {
                sigmaMatrix(i,j) = (*BIG_SIGMA)(i,j) + distribution(generator) * 0.005; // add epsilon to SIGMA
            }
        }
        */

	//int nT = omp_get_num_procs();
	//omp_set_num_threads(1);
	//if (!lowrank) {
	//	printf("safe to invert\n");
        invSigmaMatrix = inv(sigmaMatrix);
        sigmaDet       = det(sigmaMatrix);
	//} else {
	//	printf("don't invert\n");
        //  invSigmaMatrix = nullptr;
	//  sigmaDet = 0;
	//}
	//omp_set_num_threads(nT);

    }

    ~PostCal() {
        delete [] histValues;
        delete [] postValues;
	delete [] noCausal;
	delete [] sharedPips;
	delete [] sharedLL;
	delete [] notSharedLL;
    }

    // https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933#72073933
// not sure about the noexcept
struct VecHash {
    std::size_t operator()(std::vector<int> const& vec) const noexcept {
        std::size_t seed = vec.size();

        for (int v : vec) {
            uint32_t x = static_cast<uint32_t>(v);

            x = ((x >> 16) ^ x) * 0x45d9f3bu;
            x = ((x >> 16) ^ x) * 0x45d9f3bu;
            x = (x >> 16) ^ x;
            seed ^= std::size_t(x) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        }

        return seed;
    }
};

    /*
     * Calculate prior probability of given configuration vector
     * */
    double log_prior(vector<int> configure, int numCaual, int **causal_bool_per_study);

    /*
     construct sigma_C by the kronecker product in paper, it is mn by mn. the variance for vec(lambdaC)|vec(C)
     :param configure the causal status vector of 0 and 1
     :return diagC is the variance matrix for (lamdaC|C)
     */
    mat construct_diagC(vector<int> configure, int numCausal, int **causal_idx_per_study, int **causal_bool_per_study);

    /*
     compute likelihood of each configuration by Woodbury
     :param configure the causal status vector of 0 and 1
     :param stat the z-score of each snp
     :param sigma_g_squared the non-centrality param
     :return likelihood of the configuration
     */
    double likelihood(vector<int> configure, vector<double> * stat, double sigma_g_squared, mat sigmaC) ;

    /*
     compute likelihood for low rank matrices
     :param configure the causal status vector of 0 and 1
     :param stat the z-score of each snp
     :param sigma_g_squared the non-centrality param
     :return likelihood of the configuration
     */
    double lowrank_likelihood(vector<int> configure, vector<double> * stat, double sigma_g_squared, mat sigmaC) ;

    /*
     find the next binary configuration based on the previous config and size of vector
     */
    int nextBinary(vector<int>& data, int size) ;

    vector<vector<int>> get_nbdplus(vector<int> curr_config);
    vector<vector<int>> get_nbdminus(vector<int> curr_config);
    vector<vector<int>> get_nbdzero(vector<int> curr_config)

    /*
     find the total likelihood given the z_score and sigma_g_squared
     */
    double computeTotalLikelihood(vector<double> * stat, double sigma_g_squared) ;
    double computeTotalLikelihoodGivenConfigs(vector<double> * stat, double sigma_g_squared) ;

    bool checkOR(int **causal_bool_per_study_for_config, const int num_of_studies, int numCausal);
    bool checkAND(int **causal_bool_per_study_for_config, const int num_of_studies, int numCausal);

    /*find configuration from iteration in string*/
    vector<int> findConfig(int iter);

    vector<int> constructConfig(vector<int> input_causal_locs);

    /*
     greedy algorithm to find minimal set
     @param stat is the z-scores
     @param sigma is the correaltion matrix
     */
    vector<char> findOptimalSetGreedy(vector<double> * stat, double sigma_g_squared, vector<int> *rank,  double inputRho, string outputFileName, double cutoff_threshold);

    /*
     print the hist file, which is the likelihood of the set containing 0, 1, 2... up to the number of max snps
     */
    void printHist2File(string fileName) {
        exportVector2File(fileName, histValues, maxCausalSNP+1);
    }
    
    double special_exp(double post, double total) {
      if ( post == 0 ) {
        return 0;
      } else {
        return exp(post - total);
      }
    }

    /*
     print to the .post file as well as no causal
     */
    void printPost2File(string fileName) {
        double total_post = 0;
        for(int i = 0; i < totalSnpCount; i++) {
            total_post = addlogSpace(total_post, postValues[i]);
	}
	total_post = totalLikeLihoodLOG;


        int start_offset = 0;
	int end_offset = num_snps_all[0];
	for ( int s = 0; s < num_of_studies; s++ ) {
          ofstream outputFile;
          string outFileNameSet = string(fileName)+"_study"+std::to_string(s)+"_post.txt";
          outputFile.open(outFileNameSet.c_str());
          outputFile << "SNP_ID\tProb_in_pCausalSet" << endl;
          int j = 0;
          for(int i = start_offset; i < end_offset; i++) {
            outputFile << (*SNP_NAME)[s][j] << "\t" << special_exp(postValues[i], total_post)  << endl; 
	    j += 1;
          }
	  start_offset = end_offset;
	  if ( s!= num_of_studies - 1 ) {
	    end_offset += num_snps_all[s+1];
	  }
	  outputFile.close();
	}
	ofstream outputFile;
	string outFileNameSet = string(fileName)+"_nocausal.txt";
        outputFile.open(outFileNameSet.c_str());
        for ( int s = 0; s < num_of_studies; s++ ) {
          outputFile << special_exp(noCausal[s], total_post) << endl;
        }
        outputFile.close();

	string outFileNameShared = string(fileName)+"_shared_pips.txt";
	outputFile.open(outFileNameShared.c_str());
	outputFile << "SNP_ID\tshared_pip\tshared_ll\tnotshared_ll" << endl;
	for ( int i = 0; i < unionSnpCount; i++ ) {
          outputFile << all_snp_pos[i] << "\t" << special_exp(sharedPips[i],total_post) << "\t" << sharedLL[i] << "\t" << notSharedLL[i] << endl;
	}
	outputFile.close();

        /* old post file output
        outfile << "SNP_ID\tProb_in_pCausalSet\tCausal_Post._Prob." << endl;
        for(int i = 0; i < snpCount; i++) {
            outfile << (*SNP_NAME)[0][i] << "\t" << exp(postValues[i]-total_post) << "\t" << exp(postValues[i]-totalLikeLihoodLOG) << endl;
        }
        */
    }

};

#endif
