#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <chrono>

#include <armadillo>

#include "postcal.h"

using namespace std;
using namespace arma;

class Model{
public:
    double sharing_param;
    double rho;
    double gamma;
    int snpCount;
    const int totalCausalSNP;
    vector<mat> * sigma;
    vector< vector<double> > * z_score;
    vector<char> * pcausalSet;
    vector<int> * rank;
    bool histFlag;  // to out the probaility of different number of causal SNP
    PostCal * post;
    vector< vector<string> > * snpNames;
    vector<string> ldDir;
    vector<string> zDir;
    string snpMapFile;
    string configsFile;
    int num_configs;
    int num_groups;
    vector<int> sample_sizes;
    vector<int> num_causal; //number of causal snps in each study
    vector<int> num_snps_all; //number of snps in each study
    string outputFileName;
    double tau_sqr;
    double sigma_g_squared;
    const int num_of_studies;
    vector<double> S_LONG_VEC;
    bool haslowrank = false;
    double cutoff_threshold;
    vector<vector<int>> idx_to_snp_map; //union position to study level position
    vector<vector<int>> idx_to_union_pos_map; //study level position to union position
    vector<string> all_snp_pos;

    /*
     consrtuctor for Model
     */
    Model(vector<string> ldDir, vector<string> zDir, string snpMapFile, string configsFile, int num_configs, int num_groups, vector<int> sample_sizes, vector<int> num_causal, string outputFileName, const int totalCausalSNP, double sharing_param, double rho, bool histFlag, double gamma=0.01, double tau_sqr = 0.2, double sigma_g_squared = 5.2, double cutoff_threshold = 0) : totalCausalSNP(totalCausalSNP), num_of_studies(ldDir.size()) {
        this->histFlag = histFlag;
	this->sharing_param = sharing_param;
        this->rho = rho;
        this->gamma = gamma;
        this->ldDir = ldDir;
        this->zDir  = zDir;
	this->snpMapFile = snpMapFile;
	this->configsFile = configsFile;
	this->num_configs = num_configs;
	this->num_groups = num_groups;
        this->outputFileName = outputFileName;
//        this->totalCausalSNP = totalCausalSNP;
        this->tau_sqr = tau_sqr;
        this->sigma_g_squared = sigma_g_squared;
        this->sample_sizes = sample_sizes;
	this-> num_causal = num_causal;
        this->cutoff_threshold = cutoff_threshold;

        //fileSize(ldFile, tmpSize);
        sigma      = new vector<mat>;
        z_score    = new vector<vector<double> >;
        snpNames   = new vector<vector<string> >;


        for(int i = 0; i < ldDir.size(); i++) {
            string ld_file = ldDir[i];
            string z_file = zDir[i];

            vector<double>* temp_LD = new vector<double>;
            vector<string> temp_names;
            vector<double> temp_z;

            importData(ld_file, temp_LD);
            importDataFirstColumn(z_file, temp_names);
            importDataSecondColumn(z_file, temp_z);

            int numSnps = sqrt(temp_LD->size());
	    num_snps_all.push_back(numSnps);
	    if (numSnps != temp_names.size()) {
              printf("ERROR: LD matrix is size %d x %d but zscores has %lu snps\n. Check LD file for nans.\n", numSnps, numSnps, temp_names.size());
	      exit(1);
	    }
	    printf("pushing back num snps %d for study %d\n", i, numSnps);
            mat temp_sig;
            temp_sig = mat(numSnps, numSnps);
            for (int i = 0; i < numSnps; i++){
                for (int j = 0; j< numSnps; j++){
                    temp_sig(i,j) = temp_LD->at(i * numSnps + j);
                }
            }

	    vector<int> idx_to_snp_studyi;
	    idx_to_snp_map.push_back(idx_to_snp_studyi);
	    vector<int> idx_to_union_pos_studyi;
	    idx_to_union_pos_map.push_back(idx_to_union_pos_studyi);

            sigma->push_back(temp_sig);
            snpNames->push_back(temp_names);
            z_score->push_back(temp_z);

            delete temp_LD;
        }
	printf("len of snpnames is %ld\n", snpNames->size());
	for ( int i = 0; i < 2; i++ ) {
          printf("len of snpnames vector %d is %ld\n", i, (*snpNames)[i].size());
	}


//        num_of_studies(snpNames->size());

	importSnpMap(snpMapFile, num_of_studies+1, &all_snp_pos, &idx_to_snp_map);
	
	for ( int i = 0; i < num_of_studies; i++ ) {
          for ( int j = 0; j < idx_to_snp_map[i].size(); j++ ) {
            if (idx_to_snp_map[i][j] >= 0 ) {
              idx_to_union_pos_map[i].push_back(j);
	    }
	  }
          if ( idx_to_union_pos_map[i].size() != num_snps_all[i] ) {
            printf("Invariant does not hold\n");
            exit(1);
	  }
	}

	int totalSnpCount = std::accumulate(num_snps_all.begin(), num_snps_all.end(), 0);

        snpCount = (*snpNames)[0].size();
        pcausalSet = new vector<char>(totalSnpCount,'0');


        rank = new vector<int>(totalSnpCount, 0);

        for (int i = 0; i < z_score->size(); i++){
            for(int j = 0; j < (*z_score)[i].size(); j++){
                S_LONG_VEC.push_back((*z_score)[i][j]);
            }
        }

        /* sigma_g_squared is set to max(5.2, max(abs(z-score)))
        for(int i = 0 ; i < num_of_studies; i++){
            for (int j = 0; j < snpCount; j++){
                if(abs(S_LONG_VEC.at(i*snpCount + j)) > sigma_g_squared){
                    sigma_g_squared = abs(S_LONG_VEC.at(i*snpCount + j));
                }
            }
        }
        */

        //make positive definite
        for (int i = 0; i < sigma->size(); i++){
            //check for low rank
            //if(arma::rank(sigma->at(i)) < snpCount){
	    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		    //printf("rank of study %d is %lld\n", i, arma::rank(sigma->at(i)));
		    //printf("num snps study %d is %d\n", i, num_snps_all[i]);
            //if(arma::rank(sigma->at(i)) < num_snps_all[i]){
	    int min_num_snps = num_snps_all[0];
	    for ( int num : num_snps_all ) {
	      if ( num < min_num_snps ) {
	        min_num_snps = num;
	      }
	    }
            if(arma::rank(sigma->at(i)) < min_num_snps){
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time to check rank = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
                haslowrank = true;
                std::cout << "study " << i << " has low rank. Implementing low_rank method.\n";
            }
            haslowrank = true;
            
            //makeSigmaPositiveSemiDefinite(&(sigma->at(i)), snpCount);
	    begin = std::chrono::steady_clock::now();
            makeSigmaPositiveSemiDefinite(&(sigma->at(i)), num_snps_all[i]);
	    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	    std::cout << "Time to make psd = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
        }

        //mat* BIG_SIGMA = new mat(snpCount * num_of_studies, snpCount * num_of_studies, fill::zeros);
	int num_total_snps = std::accumulate(num_snps_all.begin(), num_snps_all.end(), 0);
        mat* BIG_SIGMA = new mat(num_total_snps, num_total_snps, fill::zeros);
        for (int i = 0 ; i < num_of_studies; i++){
            //mat temp_sigma = mat(num_of_studies , num_of_studies, fill::zeros);
            //temp_sigma(i,i) = 1;
            //temp_sigma = kron(temp_sigma, sigma->at(i));
            //(*BIG_SIGMA) = (*BIG_SIGMA) + temp_sigma;
	    int sum_msubj_until_i = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
	    (*BIG_SIGMA).submat(sum_msubj_until_i, sum_msubj_until_i, sum_msubj_until_i+num_snps_all[i]-1, sum_msubj_until_i+num_snps_all[i]-1) = sigma->at(i);
        }
        
        //if low rank, BIG_SIGMA = BIG_B, Stat matrix has new distribution
	//omp_set_num_threads(1);
        if(haslowrank == true){
            //construct big B
            mat* BIG_B = new mat(num_total_snps, num_total_snps, fill::zeros);
            for(int i = 0; i<num_of_studies; i++){
                mat* tmpmat = new mat(num_snps_all[i], num_snps_all[i], fill::zeros);
		int sum_msubj_until_i = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
                //*tmpmat = BIG_SIGMA->submat(i*snpCount,i*snpCount,(i+1)*snpCount-1,(i+1)*snpCount-1);                
                *tmpmat = BIG_SIGMA->submat(sum_msubj_until_i, sum_msubj_until_i, sum_msubj_until_i+num_snps_all[i]-1, sum_msubj_until_i+num_snps_all[i]-1);
                mat* tmpOmega = new mat(num_snps_all[i],num_snps_all[i],fill::zeros);
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                tmpOmega = eigen_decomp(tmpmat,num_snps_all[i]);
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time for eigen decomp = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

                *tmpOmega = abs(*tmpOmega);

                //construct B
                mat trans_Q = trans(*tmpmat);
                mat sqrt_Omega = sqrt(*tmpOmega);
                mat B_each = sqrt_Omega * trans_Q;

                //merge to Big B
                //mat temp_b = mat(num_of_studies , num_of_studies, fill::zeros);
                //temp_b(i,i) = 1;
                //temp_b = kron(temp_b, B_each);
                //(*BIG_B) = (*BIG_B) + temp_b;
		(*BIG_B).submat(sum_msubj_until_i, sum_msubj_until_i, sum_msubj_until_i+num_snps_all[i]-1, sum_msubj_until_i+num_snps_all[i]-1) = B_each;

                //update S_LONG_VEC
                mat* z_score = new mat(num_snps_all[i],1,fill::zeros);
                for(int j = 0; j < num_snps_all[i]; j++){ 
                    (*z_score)(j,0) = S_LONG_VEC[sum_msubj_until_i+j]; //msubj the j does not correspond to j in this for loop, it is just an idx
                }

		//int nT = omp_get_num_procs();
		//omp_set_num_threads(1);
                mat tmpS = inv(sqrt_Omega) * trans_Q;
		//omp_set_num_threads(nT);
                mat lowS = tmpS * (*z_score);

                for(int j = 0; j < num_snps_all[i]; j++){
                    S_LONG_VEC[sum_msubj_until_i+j] = lowS(j,0);
                }
                delete(tmpmat);
                delete(tmpOmega);
                delete(z_score);
            }

	    //delete(BIG_SIGMA);
            *BIG_SIGMA = *BIG_B;
            delete(BIG_B);
        }
        post = new PostCal(BIG_SIGMA, &S_LONG_VEC, snpCount, configsFile, num_configs, num_groups, totalCausalSNP, num_causal, snpNames, sharing_param, gamma, tau_sqr, sigma_g_squared, num_of_studies, sample_sizes, num_snps_all, haslowrank, idx_to_snp_map, idx_to_union_pos_map, all_snp_pos);
    }

    /*
     run the greedy algorithm
     @param no param
     @return no return
     */
    void run() {
        (*pcausalSet) = post->findOptimalSetGreedy(&S_LONG_VEC, sigma_g_squared, rank, rho, outputFileName, cutoff_threshold);
    }

    /*
     finish by by printing the set, post and hist file
     @param no param
     @return no return
     */
    void finishUp() {
	//printf("causal set:\n");
	//printCharVec(*pcausalSet);
	int start_offset = 0;
	int end_offset = num_snps_all[0];
	for ( int s = 0; s < num_of_studies; s++ ) {
          ofstream outputFile;
          string outFileNameSet = string(outputFileName)+"_study"+std::to_string(s)+"_set.txt";
          outputFile.open(outFileNameSet.c_str());
          int j = 0;
          for(int i = start_offset; i < end_offset; i++) {
            if((*pcausalSet)[i] == '1') {
                outputFile << (*snpNames)[s][j] << endl;
	    }
	    j += 1;
          }
	  start_offset = end_offset;
	  if ( s != num_of_studies-1 ) {
	    end_offset += num_snps_all[s+1];
	  }
	  outputFile.close();
	}
        post->printPost2File(string(outputFileName));
	/* commenting out for now
        //outputs the histogram data to file
        if(histFlag)
            post->printHist2File(string(outputFileName)+"_hist.txt");
	    */
    }

    // destructor
    ~Model() {
        delete z_score;
        delete sigma;
        delete snpNames;
        delete pcausalSet;
        delete rank;
	//TODO need to delete BIG SIGMA
    }
};

#endif
