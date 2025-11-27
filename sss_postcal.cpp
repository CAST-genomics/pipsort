#include <vector>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <iostream>
#include <armadillo>
#include <iomanip>
#include <vector>
#include <math.h>
#include "util.h"
#include "postcal.h"
#include <sys/mman.h>
#include <omp.h>
#include <ctime>
#include <chrono>

using namespace arma;

double PostCal::log_prior(vector<int> configure, int numCausal, int **causal_bool_per_study) {
  if (num_of_studies > 2) {
    cout << "This prior does not work for more than 2 studies yet\n";
    exit(1);
  } 
  double p_of_c = 0;

  
  if ( sharing_param != 0 ) {
    for ( int i = 0; i < numCausal; i++ ) {
      if ( causal_bool_per_study[0][i] == causal_bool_per_study[1][i] ) {
        if ( causal_bool_per_study[0][i] == 1 ) {
          p_of_c += log(sharing_param);
 	  //printf("shared causal\n");
          } else {
          cout << "This case in prior should not happen\n"; //this is the all zeros case
	  exit(1);
        }
      } else {
        p_of_c += log((1-sharing_param)*0.5); //0.5 adjustment for 2 studies bc not shared causal has two possibilities
	//printf("not shared causal\n");
      }
    }
  }
  //printf("partial prior is %f\n", p_of_c);

  int num_1_in_either = 0;
  for ( int i = 0; i < numCausal; i++ ) {
    if (causal_bool_per_study[0][i] == 1 or causal_bool_per_study[1][i] == 1 ) {
      p_of_c += log(gamma);
      num_1_in_either += 1;
      //printf("causal snp\n");
    } 
  }
  p_of_c += ( unionSnpCount - num_1_in_either ) * log(1 - gamma);
  //printf("num not causal snp = %d\n", unionSnpCount - num_1_in_either);
  //printf("prior is %f\n", p_of_c);
  //printf("num of 1 in either is %d\n", num_1_in_either);
  
  return p_of_c;
}


// calibrate for sample size imbalance
mat PostCal::construct_diagC(vector<int> configure, int numCausal, int **causal_idx_per_study, int **causal_bool_per_study) {
    mat Identity_M = mat(num_of_studies, num_of_studies, fill::eye);
    mat Matrix_of_sigmaG = mat(num_of_studies, num_of_studies);
    int min_size = * std::min_element(sample_sizes.begin(), sample_sizes.end());

    /*for (int i = 0; i < num_of_studies; i ++) {
        for (int j = 0; j < num_of_studies; j ++) {
            if (i == j) // diagonal: scaled variance
                Matrix_of_sigmaG(i, j) = s_squared * (double(sample_sizes[i]) / min_size);
            else // off-diagonal: covariance
                Matrix_of_sigmaG(i, j) = s_squared * sqrt(long(sample_sizes[i]) * long(sample_sizes[j])) / min_size;
        }
    } */   
    //mat temp1 = t_squared * Identity_M + Matrix_of_sigmaG;
    //mat temp1 = Matrix_of_sigmaG;
    //mat temp2 = mat(snpCount, snpCount, fill::zeros);
    //for(int i = 0; i < snpCount; i++) {
    //    if (configure[i] == 1)
    //       temp2(i, i) = 1;
    //}
    vec diagC_main_diag(totalSnpCount, fill::zeros);
    for ( int i = 0; i < num_of_studies; i++ ) {
      for ( int j = 0; j < numCausal; j++ ) {
        if ( causal_bool_per_study[i][j] == 1 ) {
          int offset_studyi = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
	  int loc_in_studyi = causal_idx_per_study[i][j];
	  diagC_main_diag[offset_studyi + loc_in_studyi] = s_squared * (double(sample_sizes[i]) / min_size) + t_squared;
	}
      }
    }

    mat diagC = diagmat(diagC_main_diag);
    //printf("diagonal diagC\n");
    //diagC.print(std::cout);
    //mat diagC = kron(temp1, temp2);
    
    //adjust off diagonals to have sigma^2 when commonly causal in two studies
    for ( int s1 = 0; s1 < num_of_studies; s1++ ) {
      for ( int s2 = 0; s2 < num_of_studies; s2++ ) {
        for ( int j = 0; j < numCausal; j++ ) {
            if ( s1 == s2 ) { continue; } //already took care of diagonals
	    if (causal_bool_per_study[s1][j] == causal_bool_per_study[s2][j] && causal_bool_per_study[s1][j] == 1) { //causal in both
	      int study1_offset = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+s1, 0);    
	      int study2_offset = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+s2, 0);    
	      int loc_in_study1 = causal_idx_per_study[s1][j];
	      int loc_in_study2 = causal_idx_per_study[s2][j];
	      int s1_idx = study1_offset + loc_in_study1;
	      int s2_idx = study2_offset + loc_in_study2;
	      diagC(s1_idx, s2_idx) = s_squared * sqrt(long(sample_sizes[s1]) * long(sample_sizes[s2])) / min_size;
	      diagC(s2_idx, s1_idx) = s_squared * sqrt(long(sample_sizes[s1]) * long(sample_sizes[s2])) / min_size;
	  }
	}
      }
    }
    
    //printf("diagC again\n");
    //diagC.print(std::cout);
    return diagC;
}

double PostCal::likelihood(vector<int> configure, vector<double> * stat, double sigma_g_squared, mat sigmaC) {
    //int causalCount = 0;
    int totalCausalCount = 0;
    double matDet   = 0;
    double res      = 0;

    for(int i = 0; i < totalSnpCount; i++) {
        totalCausalCount += configure[i];
    }
    if (totalCausalCount == 0) {
      cout << "This should have been taken care of in total likelihood function\n";
      exit(1);
    }
    /*if(totalCausalCount == 0){ //TODO how to handle this case for multiple causal vectors
        mat tmpResultMatrixNM = statMatrixtTran * invSigmaMatrix;
        mat tmpResultMatrixNN = tmpResultMatrixNM * statMatrix;

        res = tmpResultMatrixNN(0,0);
        matDet = sigmaDet;
        return (-res/2-sqrt(abs(matDet)));
    }*/

    //mat sigmaC = construct_diagC(configure);
    int index_C = 0;
    mat sigmaMatrixTran = sigmaMatrix.t();

    // U is kn by mn matrix of columns corresponding to causal SNP in sigmacC
    // In unequal sample size studies, U is adjusted for the sample sizes
    //mat U(causalCount * num_of_studies, snpCount * num_of_studies, fill::zeros);
    mat U(totalCausalCount, totalSnpCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                U(index_C, j) = sigmaC(i, j);
            }
            index_C ++;
        }
    }
    
    index_C = 0;
    
    // V is mn by kn matrix of rows corresponding to causal SNP in sigma
    // In unequal sample size studies, V does not change
    //mat V(causalCount * num_of_studies, snpCount * num_of_studies, fill::zeros);
    mat V(totalCausalCount, totalSnpCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                V(index_C, j) = sigmaMatrixTran(i, j);
            }
            index_C ++;
        }
    }
    V = V.t();

    // UV = SigmaC * Sigma (kn by kn)
    mat UV(totalCausalCount, totalCausalCount, fill::zeros);
    UV = U * V;

    //mat I_AA   = mat(snpCount, snpCount, fill::eye); //can ignore, does not get used here
    mat tmp_CC = mat(totalCausalCount, totalCausalCount, fill::eye) + UV;
    matDet = det(tmp_CC) * sigmaDet;

    mat temp1 = invSigmaMatrix * V;
    mat temp2 = mat(totalSnpCount, totalCausalCount, fill::zeros);
    
    #pragma omp critical
    temp2 = temp1 * pinv(tmp_CC);

    mat tmp_AA = invSigmaMatrix - temp2 * U ;

    mat tmpResultMatrix1N = statMatrixtTran * tmp_AA;
    mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
    res = tmpResultMatrix11(0,0);

    if(matDet==0) {
        cout << "Error the matrix is singular and we fail to fix it. (reg lkl)" << endl;
        exit(0);
    }

    /*
     We compute the log of -res/2-log(det) to see if it is too big or not.
     In the case it is too big we just make it a MAX value.
     */
    double tmplogDet = log(sqrt(abs(matDet)));
    double tmpFinalRes = -res/2 - tmplogDet;

    return tmpFinalRes;
}

//here we still use Woodbury matrix, here sigma_matrix is B, and S is updated already
double PostCal::lowrank_likelihood(vector<int> configure, vector<double> * stat, double sigma_g_squared, mat sigmaC) {
    //int causalCount = 0;
    int totalCausalCount = 0;
    double matDet   = 0;
    double res      = 0;

    for(int i = 0; i < totalSnpCount; i++) {
        totalCausalCount += configure[i];
    }
    if (totalCausalCount == 0) {
      cout << "This should have been taken care of in total likelihood function\n";
      exit(1);
    }
    
    /*if(totalCausalCount == 0){ //TODO ok for now, but need to think
        mat tmpResultMatrixNN = statMatrixtTran * statMatrix;
        res = tmpResultMatrixNN(0,0);
        matDet = 1;
        return (-res/2-sqrt(abs(matDet)));
    }*/

    //mat sigmaC = construct_diagC(configure);
    int index_C = 0;
    mat sigmaMatrixTran = sigmaMatrix.t();

    // In unequal sample size studies, U is adjusted for the sample sizes
    // here we make U = B * sigmaC, this is still kn by mn
    mat U = mat(totalCausalCount, totalSnpCount, fill::zeros);

    mat small_sigma(totalSnpCount, totalCausalCount, fill::zeros);
    mat small_sigmaC(totalCausalCount, totalCausalCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                small_sigma(j, index_C) = sigmaMatrix(j, i);
            }
            small_sigmaC(index_C, index_C) = sigmaC(i, i);
            index_C++;
        }
    }
    U = small_sigma * small_sigmaC;
    U = U.t();

    index_C = 0;

    // here V is B_trans, this is mn by kn
    mat V = mat(totalCausalCount, totalSnpCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                V(index_C, j) = sigmaMatrixTran(i, j);
            }
            index_C ++;
        }
    }
    V = V.t();

    // UV = B * SigmaC * Btrans (kn by kn)
    mat UV(totalCausalCount, totalCausalCount, fill::zeros);
    UV = U * V;

    mat* I_AA   = new mat(totalSnpCount, totalSnpCount, fill::eye);
    mat tmp_CC = mat(totalCausalCount, totalCausalCount, fill::eye) + UV;
    matDet = det(tmp_CC);

    mat temp2 = mat(totalSnpCount, totalCausalCount, fill::zeros);
    #pragma omp critical
    temp2 = V * pinv(tmp_CC);

    mat tmp_AA = *I_AA - temp2 * U ;

    mat tmpResultMatrix1N = statMatrixtTran * tmp_AA;
    mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
    res = tmpResultMatrix11(0,0);

    delete(I_AA);
    
    if(matDet==0) {
        cout << "Error the matrix is singular and we fail to fix it (low rank lkl)." << endl;
        exit(0);
    }

    /*
     We compute the log of -res/2-log(det) to see if it is too big or not.
     In the case it is too big we just make it a MAX value.
     */
    double tmplogDet = log(sqrt(abs(matDet)));
    double tmpFinalRes = -res/2 - tmplogDet;

    return tmpFinalRes;
}

vector<vector<int>> PostCal::get_nbdplus(vector<int> curr_config) {
    //each vector will be the list of indices corresponding to causal snps
    vector<vector<int>> plusone_configs;
    int curr_num_causal = 0;
    vector<int> curr_causal_locs;
    for ( int i = 0; i < curr_config.size(); i++ ) {
        if curr_config[i] == 1 {
           curr_num_causal += 1
           curr_causal_locs.push_back(i)
	}
    }
    int num_plus_configs = curr_config.size() - curr_num_causal;
    for ( int i = 0; i < curr_config.size(); i++ ) {
      if curr_config[i] == 0 {
          vector<int> new_config;
	  new_config.push_back(curr_config[i]);
	  new_config.insert(new_config.end(), curr_causal_locs.begin(), curr_causal_locs.end());
	  plusone_configs.push_back(new_config);
      }
    }
    if ( plusone_configs.size() != num_plus_configs ) {
       printf("plusone configs WRONG size\n");
    }
    return plusone_configs;
}

vector<vector<int>> PostCal::get_nbdminus(vector<int> curr_config) {
    //each vector will be the list of indices corresponding to causal snps
    vector<vector<int>> minusone_configs;
    int curr_num_causal = 0;
    vector<int> curr_causal_locs;
    for ( int i = 0; i < curr_config.size(); i++ ) {
        if curr_config[i] == 1 {
           curr_num_causal += 1
           curr_causal_locs.push_back(i)
	}
    }
    int num_minus_configs = curr_num_causal;
    for ( int i = 0; i < num_minus_configs; i++ ) {
        vector<int> new_config;
	for ( int j = 0; j < curr_num_causal; j++ ) {
            if ( i != j ) {
                new_config.push_back(curr_causal_locs[j]);
	    }
	}
        minusone_configs.push_back(new_config);
    }
    if ( minusone_configs.size() != num_minus_configs ) {
       printf("minusone configs WRONG size\n");
    }
    return minusone_configs;
}


vector<vector<int>> PostCal::get_nbdzero(vector<int> curr_config) {
    //each vector will be the list of indices corresponding to causal snps
    vector<vector<int>> nbdzero_configs;
    int curr_num_causal = 0;
    vector<int> curr_causal_locs;
    for ( int i = 0; i < curr_config.size(); i++ ) {
        if curr_config[i] == 1 {
           curr_num_causal += 1
           curr_causal_locs.push_back(i)
	}
    }
    int num_nbdzero_configs = curr_num_causal * (curr_config.size() - curr_num_causal);
    vector<vector<int>> minusone_configs = get_nbdminus(curr_config); 
    for ( int i = 0; i < curr_config.size(); i++ ) {
	if ( curr_config[i] == 0 ) {
	    for ( vector<int> v : minusone_configs ) {
		vector<int> new_config;
		new_config.pushback(i);
		new_config.insert(new_config.end(), v.begin(), v.end());
		nbdzero_configs.pushback(new_config);
	    }
	}
    }
    if ( nbdzero_configs.size() != num_nbdzero_configs ) {
       printf("nbdzero configs WRONG size\n");
    }
    return nbdzero_configs;
}


double PostCal::sss_computeTotalLikelihood(vector<double>* stat, double sigma_g_squared) {
    double sumLikelihood = 0;
    int mycount = 0;
    printf("num total configs = %d\n", mycount);

    int unionSnpCount = all_snp_pos.size();
    cout << "Max Causal = " << maxCausalSNP << endl;
    cout << "Union Snp Count = " << unionSnpCount << endl;

    //store likeliehood of explored configs
    std::unordered_map<vector<int>, double> config_hashmap;

    //clock_t start = clock();
    //TODO initialized as this for now
    vector<int> configure(unionSnpCount, 0);
    int num;

    int curr_iter = 0;


/*
 some rough pseudocode
 init configure somehow
 will need a hash map to keep track of configs I've explored
 key for hash map can be indices in confif (size=unionsnpcount) that are causal
 for t = 1, ... T
    get nbd(configure) = gammaplus, gamma0, gammaminus
    for each c in nbd(configure):
       expand into totalConfig and eval likelihood, keeping track of max likelihood
       add c and it's max likelihood to hash map
    then apply sampling scheme?

 *
 *
 */


    
    int nP = omp_get_num_procs();
    printf("num procs = %d\n", nP);
    printf("max num procs = %d\n", omp_get_max_threads());
    printf("num threads = %d\n", omp_get_num_threads());
    int*** thread_mem_idx = (int***)malloc(nP * sizeof(int**));
    int*** thread_mem_bool = (int***)malloc(nP * sizeof(int**));
    int*** thread_mem_bool_for_config = (int***)malloc(nP * sizeof(int**));
    for ( int b = 0; b < nP; b++ ) {
      thread_mem_idx[b] = (int**)malloc(num_of_studies * sizeof(int*));
      thread_mem_bool[b] = (int**)malloc(num_of_studies * sizeof(int*));
      thread_mem_bool_for_config[b] = (int**)malloc(num_of_studies * sizeof(int*));
      for ( int q = 0; q < num_of_studies; q++ ) {
        thread_mem_idx[b][q] = (int*)malloc(3 * sizeof(int));
        thread_mem_bool[b][q] = (int*)malloc(3 * sizeof(int));
        thread_mem_bool_for_config[b][q] = (int*)malloc(3 * sizeof(int));
      }
    }
    
    omp_set_num_threads(nP);

    //TODO hash for vector<int>
    //TODO move the pragma somewhere else
    //#pragma omp parallel for schedule(static,chunksize) private(configure,num)
    int total_iteration = 10; //TODO for now
    for(int iter = 0; iter < total_iteration; iter++) {

	int numCausal = 0;
	vector<int> causal_locs;
    	for ( int i = 0; i < configure.size(); i++ ) {
            if configure[i] == 1 {
               numCausal += 1
               causal_locs.push_back(i)
            }
        }
	if ( numCausal == 0 ) {
            causal_locs.push_back(-1); //repr for all zeros config
	}

	vector<vector<int>> nbdzero = get_nbdzero(configure);
	vector<vector<int>> nbdminus = get_nbdminus(configure);
	vector<vector<int>> nbdplus = get_nbdplus(configure);
	int num_zero = nbdzero.size();
	int num_minus = nbdminus.size();
	int num_plus = nbdplus.size();

	vector<vector<int>> nbd;
        nbd.reserve(nbdzero.size() + nbdminus.size() + nbdplus.size());

        nbd.insert(nbd.end(),
            std::make_move_iterator(nbdzero.begin()),
            std::make_move_iterator(nbdzero.end()));

        nbd.insert(nbd.end(),
            std::make_move_iterator(nbdminus.begin()),
            std::make_move_iterator(nbdminus.end()));

        nbd.insert(nbd.end(),
            std::make_move_iterator(nbdplus.begin()),
            std::make_move_iterator(nbdplus.end()));



	bool make_updates = true;	
	//check if we have already computed likeliehood for this config
	if (config_hashmap.find(causal_locs) != config_hashmap.end()) {
	    make_updates = false;
        }
	
	double tmp_likelihood = expand_and_compute_lkl(configure, make_updates);
	vector<double> probs(nbd.size(), 0);
	//TODO pragma here
	//change from for each to just for, change push_back to probs[i] = X
	for ( int i = 0; i < nbd.size(); i++ ) {
	   double v_lkl = expand_and_compute_lkl(v, false); //compute lkl but no update
	   probs[i] = v_lkl;
	}

	//sampling

	double weight_zero = 0.0; size_t zero_sample = nbd.size(); //nbd.size() is an invalid index
	double weight_minus = 0.0; size_t minus_sample = nbd.size();
	double weight_plus = 0.0; size_t plus_sample = nbd.size();

	std::random_device rd;
	std::mt19937 gen(rd());
	if ( num_zero != 0 ) {
	    std::discrete_distribution<size_t> dist(probs.begin(), probs.begin()+num_zero);
	    zero_sample = dist(gen);
	    weight_zero = std::accumulate(probs.begin(), probs.begin()+num_zero, 0.0);
	}
	if ( num_minus != 0 ) {
	    std::discrete_distribution<size_t> dist(probs.begin()+num_zero, probs.begin()+num_zero+num_minus);
	    minus_sample = dist(gen);
	    weight_minus = std::accumulate(probs.begin()+num_zero, probs.begin()+num_zero+num_minus, 0.0);
	}
	if ( num_zero != 0 ) {
	    std::discrete_distribution<size_t> dist(probs.begin()+num_zero+num_minus, probs.end());
	    plus_sample = dist(gen);
	    weight_plus = std::accumulate(probs.begin()+num_zero+num_minus, probs.end(), 0.0);
	}

	std::discrete_distribution<size_t> dist({weight_zero, weight_minus, weight_plus});
	size_t final_sample = dist(gen);
	configure = nbd[final_sample];

    }

    omp_set_num_threads(1);

    
    for ( int i = 0; i < nP; i++ ) {
      for ( int j = 0; j < num_of_studies; j++ ) {
        free(thread_mem_idx[i][j]);
        free(thread_mem_bool[i][j]);
        free(thread_mem_bool_for_config[i][j]);
      }
      free(thread_mem_idx[i]);
      free(thread_mem_bool[i]);
      free(thread_mem_bool_for_config[i]);
    }
    free(thread_mem_idx);
    free(thread_mem_bool);
    free(thread_mem_bool_for_config);
    

    //cout << "\ncomputing likelihood of all configurations took  " << (float)(clock()-start)/CLOCKS_PER_SEC << "seconds.\n";

    for(int i = 0; i <= maxCausalSNP; i++) { //TODO what is this for, do I need to change it
        histValues[i] = exp(histValues[i]-sumLikelihood);
    }
    printf("num total configs = %d\n", mycount);


    return(sumLikelihood);
}

double PostCal::expand_and_compute_lkl(vector<int> configure, bool make_updates) {

	int numCausal = 0;
        vector<int> causal_locs;
        for ( int i = 0; i < configure.size(); i++ ) {
            if configure[i] == 1 {
               numCausal += 1
               causal_locs.push_back(i)
            }
        }

        if ( numCausal == 0 ) { //if no causal, just update sum likelihood, nothing else should change
                vector<int> tempConfigure(totalSnpCount, 0);
                double tmp_likelihood = 0;

                if(haslowrank==true){
		//mycount += 1;
		  mat tmpResultMatrixNN = statMatrixtTran * statMatrix;
                  double res = tmpResultMatrixNN(0,0);
                  double matDet = 1;
                  double lrl = (-res/2-sqrt(abs(matDet))); //lrl = low rank likelihood
                  tmp_likelihood = lrl + unionSnpCount * log(1-gamma);
                }
                else{
		  //mycount += 1;
		  mat tmpResultMatrixNM = statMatrixtTran * invSigmaMatrix;
                  mat tmpResultMatrixNN = tmpResultMatrixNM * statMatrix;

                  double res = tmpResultMatrixNN(0,0);
                  double matDet = sigmaDet;
                  double ll = (-res/2-sqrt(abs(matDet)));
                  tmp_likelihood = ll + unionSnpCount * log(1-gamma);
                }  
		if ( make_updates ) {
		   for ( int ss = 0; ss < num_of_studies; ss++ ) {
                      noCausal[ss] = addlogSpace(noCausal[ss], tmp_likelihood);
		   }
		}

                //#pragma omp critical
		//TODO where to put the pragma
		if ( make_updates ) { 
                   sumLikelihood = addlogSpace(sumLikelihood, tmp_likelihood);
		   //emplace won't make a copy of causal_locs
		   config_hashmap.emplace(causal_locs, tmp_likelihood);
		}
		return tmp_likelihood;
	}


	
	//printf("total snp count is: %d\n", totalSnpCount);
        vector<int> startConfigure(totalSnpCount, 0);

	int pid = omp_get_thread_num();
	if ((pid < 0) || (pid >= nP)) {
          cout << "Invalid PID\n" << endl;
          exit(1); // terminate with error
        }


	int** causal_idx_per_study = thread_mem_idx[pid];
	int** causal_bool_per_study = thread_mem_bool[pid];
	for ( int i = 0; i < num_of_studies; i++ ) {
	  memset(causal_idx_per_study[i], 0, 3 * sizeof(int));
	  memset(causal_bool_per_study[i], 0, 3 * sizeof(int));
	  //memset(causal_bool_per_study_for_config[i], 0, 3 * sizeof(int));
	}

	for ( int i = 0; i < num_of_studies; i++ ) {
            for ( int j = 0; j < numCausal; j++ ) {
		int loc_in_studyi = idx_to_snp_map[i][causal_locs[j]];
		//printf("loc in study %d is %d\n", i, loc_in_studyi);
		int studyi_offset = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
                if (loc_in_studyi >= 0) { //means it exists in study i
		    causal_idx_per_study[i][j] = loc_in_studyi;
		    causal_bool_per_study[i][j] = 1;
		    startConfigure[studyi_offset + loc_in_studyi] = 1;
		} else {
                    causal_idx_per_study[i][j] = -1;
		    causal_bool_per_study[i][j] = 0;
		}
	    }	    
	}
       /* 
	printf("causal index per study\n");
        for ( int i = 0; i < 2; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
		printf("%d ", causal_idx_per_study[i][j]);
	}
	printf("\n");
	}
	printf("causal bool per study\n");
        for ( int i = 0; i < 2; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
		printf("%d ", causal_bool_per_study[i][j]);
	}
	printf("\n");
	}*/

	int total_num_additional_configs = 0;
        for ( int i = 0; i < numCausal; i++ ) {
	    for ( int j = 0; j < num_of_studies; j++ ) {
                if ( causal_bool_per_study[j][i] == 1 ) {
                    total_num_additional_configs += 1;
		}
	    }
	}	
	total_num_additional_configs = (int)(pow(2, total_num_additional_configs) - 1);
//	printf("num additional = %d\n", total_num_additional_configs);


	for ( int i = 0; i < total_num_additional_configs; i++ ) {
	  //int pid = omp_get_thread_num();
	  //if ( ( pid < 0 ) || ( pid >= nP ) ) { printf("BAD PID\n"); }
	  int **causal_bool_per_study_for_config = thread_mem_bool_for_config[pid];
	  for ( int j = 0; j < num_of_studies; j++ ) {
	    memset(causal_bool_per_study_for_config[j], 0, 3 * sizeof(int));
	  }

	  int bmask = i+1;
	  vector<int> nextConfigure(totalSnpCount, 0); 
          for ( int j = 0; j < numCausal; j++ ) {
            for ( int k = 0; k < num_of_studies; k++ ) {
	      int studyk_offset = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+k, 0);
	      if ( causal_bool_per_study[k][j] == 1 ) { //if causal snp j exists in study k
	        int loc_in_studyk = causal_idx_per_study[k][j];
		int causal_val = bmask & 0x1;
		causal_bool_per_study_for_config[k][j] = causal_val;
		nextConfigure[studyk_offset + loc_in_studyk] = causal_val;
		bmask = bmask >> 1;
	      } else {
	          continue;
	      }
	    }
	  }
	  
	  /*
	printf("causal bool per study for config\n");
        for ( int ii = 0; ii < 2; ii++ ) {
        for ( int jj = 0; jj < 3; jj++ ) {
		printf("%d ", causal_bool_per_study_for_config[ii][jj]);
	}
	printf("\n");
	}*/
	  if (!checkOR(causal_bool_per_study_for_config, num_of_studies, numCausal)) {
            continue;
	  }
	  //if (!checkAND(causal_bool_per_study_for_config, num_of_studies, numCausal)) {
	  //  continue;
	  //}
	  //printf("next config to eval\n");
	  //printVec(nextConfigure);
          double tmp_likelihood = 0;
	  double just_ll = 0;
	  double just_prior = 0;
          mat sigmaC = construct_diagC(nextConfigure, numCausal, causal_idx_per_study, causal_bool_per_study_for_config);
          //printf("sigma C\n");
	  //sigmaC.print(std::cout);

          if ( haslowrank == true ) {
             //mycount += 1;
             just_ll = lowrank_likelihood(nextConfigure, stat, sigma_g_squared, sigmaC); 
             just_prior = log_prior(nextConfigure, numCausal, causal_bool_per_study_for_config); 
             tmp_likelihood = just_ll + just_prior;
             }
          else {
             //mycount += 1;
             just_ll = likelihood(nextConfigure, stat, sigma_g_squared, sigmaC);
             just_prior = log_prior(nextConfigure, numCausal, causal_bool_per_study_for_config);
             tmp_likelihood = just_ll + just_prior;
             }
        
         #pragma omp critical
         sumLikelihood = addlogSpace(sumLikelihood, tmp_likelihood);
	 //printf("sumLikelihood is %f\n", sumLikelihood);
 
	 //printf("likelihood is %f\n", tmp_likelihood);
	 //printf("exp likelihood is %f\n", exp(tmp_likelihood));

         for ( int w = 0; w < num_of_studies; w++ ) {
	   bool allZero = true;
           for ( int v = 0; v < numCausal; v++ ) {
             if (causal_bool_per_study_for_config[w][v] == 1) {
	       allZero = false;
	       break;
	     }
	   }
	   if ( allZero ) {
		   //printf("no causal in study %d\n", w);
             noCausal[w] = addlogSpace(noCausal[w], tmp_likelihood);
	   }
	 }


	 for ( int g = 0; g < numCausal; g++ ) {
           bool sharedCausal = true;
	   for ( int h = 0; h < num_of_studies; h++ ) {
             if ( causal_bool_per_study_for_config[h][g] == 0 ) {
               sharedCausal = false;
	     }
	   }
	   if ( sharedCausal ) {
             //printf("shared causal\n");
             sharedPips[causal_locs[g]] =  addlogSpace(sharedPips[causal_locs[g]], tmp_likelihood);
	     sharedLL[causal_locs[g]] = addlogSpace(sharedLL[causal_locs[g]], just_ll);
	     //printf("shared pip = %f", sharedPips[causal_locs[g]]);
	   } else {
             notSharedLL[causal_locs[g]] = addlogSpace(notSharedLL[causal_locs[g]], just_ll);
	   }
	 }

         for(int f = 0; f < totalSnpCount; f++) {
            //for(int k = 0; k < num_of_studies; k++){
                #pragma omp critical
                postValues[f] = addlogSpace(postValues[f], tmp_likelihood * nextConfigure[f]);
		//if ( nextConfigure[f] == 1 ) {
                  //printf("updating index %d\n", f);
		  //printf("added in log space %f\n", tmp_likelihood);
		  //printf("post value for %d is %f\n", f, postValues[f]);
		//}
                //}
         }        
	//printf("post values\n");
	//for ( int f = 0; f < totalSnpCount; f++) {
        //  printf("%f ", postValues[f]);
	//}
	//printf("\n");
       }
	return tmp_likelihood;
}

bool PostCal::checkAND(int **causal_bool_per_study_for_config, int num_of_studies, int numCausal) {
  int sum = 0;
  for ( int i = 0; i < numCausal; i++ ) {
    int b = 0;
    for ( int j = 0; j < num_of_studies; j++ ) {
      b += causal_bool_per_study_for_config[j][i];
    } 
    if ( b == num_of_studies ) {
      sum += 1;
    }
  }
  if (sum == numCausal) {
    return true;
  }
  return false;
}

bool PostCal::checkOR(int **causal_bool_per_study_for_config, int num_of_studies, int numCausal) {
  int sum = 0;
  for ( int i = 0; i < numCausal; i++ ) {
    int b = 0;
    for ( int j = 0; j < num_of_studies; j++ ) {
      b += causal_bool_per_study_for_config[j][i];
    } 
    if ( b > 0 ) {
      sum += 1;
    }
  }
  if (sum == numCausal) {
    return true;
  }
  return false;
}

vector<char> PostCal::findOptimalSetGreedy(vector<double> * stat, double sigma_g_squared, vector<int> * rank,  double inputRho, string outputFileName, double cutoff_threshold) {
    double total_post = double(0);

    vector<char> causalSet(totalSnpCount,'0');

    if (configsFile != "") {
      totalLikeLihoodLOG = computeTotalLikelihoodGivenConfigs(stat, sigma_g_squared);
    } else {
      totalLikeLihoodLOG = computeTotalLikelihood(stat, sigma_g_squared);
    }

    export2File(outputFileName+"_log.txt", exp(totalLikeLihoodLOG)); //Output the total likelihood to the log File
    for(int i = 0; i < totalSnpCount; i++) {
        total_post = addlogSpace(total_post, postValues[i]);
    }
    printf("\nTotal Likelihood = %e SNP=%d \n", total_post, totalSnpCount);
    total_post = totalLikeLihoodLOG;
    printf("total post as total likelihood log = %f\n", total_post);


    for ( int i = 0; i < num_of_studies; i++ ) {
      printf("no causal just value %f\n", noCausal[i]);
      printf("Prob of no causal for study %d is %f\n", i, exp(noCausal[i]-total_post));
    }

    for ( int i = 0; i < totalSnpCount; i++ ) {
      double pip = special_exp(postValues[i], total_post);
      if ( pip > 0.05 ) {
        causalSet[i] = '1';
      }
      printf("post value for index %d = %f\n", i, pip);
    }
    //printf("causal set:\n");
    //printCharVec(causalSet);

    std::vector<data> items;
    std::set<int>::iterator it;
    //output the poster to files
    //printf("posteriors\n");
    for(int i = 0; i < totalSnpCount; i++) {
        //printf("%d==>%e ",i, postValues[i]/total_likelihood);
        items.push_back(data(exp(postValues[i]- total_post), i, 0));
        //printf("%f ", exp(postValues[i]- total_post));
    }
    printf("\n");
    int start_offset = 0;
    int end_offset = 0;
    for ( int s = 0; s < num_of_studies; s++ ) {
      end_offset += num_snps_all[s]; //update end offset, needs to happen here, not where start offset is updated
      printf("start offset = %d\n", start_offset);
      printf("end offset = %d\n", end_offset);
      std::sort(items.begin()+start_offset, items.begin()+end_offset, by_number());
      printf("sort complete %d\n", s);
      for(int i = 0; i < num_snps_all[s]; i++) {
        (*rank)[start_offset+i] = items[i].index1;
      } 
      start_offset = end_offset; //update start offset
    }

    /* replaced by above for loop
    std::sort(items.begin(), items.end(), by_number());
    for(int i = 0; i < snpCount; i++)
        (*rank)[i] = items[i].index1;
    */

    double threshold = cutoff_threshold;
    /*
    if(snpCount > 30){
        threshold = 1/(double)snpCount;
    }
    else{
        threshold = 0.1/(double)snpCount;
    }
    */
    cout << "threshold is " << threshold << "\n";
    
    //reset offsets 
    start_offset = 0;
    end_offset = 0;
    for ( int s = 0; s < num_of_studies; s++ ) {
      end_offset += num_snps_all[s];
      double rho = double(0);
      int index = 0;
      while(rho < inputRho){
        rho += special_exp(postValues[(*rank)[start_offset+index]], total_post);
        if(special_exp(postValues[(*rank)[start_offset+index]], total_post) > threshold){
            //causalSet[(*rank)[start_offset+index]] = '1';
	    double pip = special_exp(postValues[start_offset+index], total_post);
	    if (pip > 0.01) {
            printf("%d %f\n", start_offset+index, pip);
	    }
        }
        index++;
	if (index >= num_snps_all[s]) {
	  //exit(1);
	  break; //credible set no longer makes sense
	}
      }
      start_offset = end_offset;
    }

    /*
    do{
        rho += exp(postValues[(*rank)[index]]-total_post);
        causalSet[(*rank)[index]] = '1';
        printf("%d %e\n", (*rank)[index], rho);
        index++;
    } while( rho < inputRho);
    */
    printf("\n");
    return(causalSet);
}
