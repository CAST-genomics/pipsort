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


int PostCal::nextBinary(vector<int>& data, int size) {
    int i = 0;
    int total_one = 0;
    int index = size-1;
    int one_countinus_in_end = 0;

    while(index >= 0 && data[index] == 1) {
        index = index - 1;
        one_countinus_in_end = one_countinus_in_end + 1;
    }

    if(index >= 0) {
        while(index >= 0 && data[index] == 0) {
            index = index - 1;
        }
    }
    if(index == -1) {
        while(i <  one_countinus_in_end+1 && i < size) {
            data[i] = 1;
            i=i+1;
        }
        i = 0;
        while(i < size-one_countinus_in_end-1) {
            data[i+one_countinus_in_end+1] = 0;
            i=i+1;
        }
    }
    else if(one_countinus_in_end == 0) {
        data[index] = 0;
        data[index+1] = 1;
    }
    else {
        data[index] = 0;
        while(i < one_countinus_in_end + 1) {
            data[i+index+1] = 1;
            if(i+index+1 >= size)
                printf("ERROR3 %d\n", i+index+1);
            i=i+1;
        }
        i = 0;
        while(i < size - index - one_countinus_in_end - 2) {
            data[i+index+one_countinus_in_end+2] = 0;
            if(i+index+one_countinus_in_end+2 >= size) {
                printf("ERROR4 %d\n", i+index+one_countinus_in_end+2);
            }
            i=i+1;
        }
    }
    i = 0;
    total_one = 0;
    for(i = 0; i < size; i++)
        if(data[i] == 1)
            total_one = total_one + 1;

    return(total_one);
}

vector<int> PostCal::findConfig(int iter) {
    int numCausal = 0;
    int temp = iter;
    int sum = 0;
    int unionSnpCount = all_snp_pos.size();
    vector<int> config(unionSnpCount, 0);
    int comb = nCr(unionSnpCount,numCausal);
    while(temp > comb) {
        temp = temp - comb;
        numCausal++;
        sum = sum + comb;
        comb = nCr(unionSnpCount,numCausal);
    }

    int times = iter - sum; //this is the number of times we use find_next_binary
    for(int i = 0; i < numCausal; i++){
        config[i] = 1;
    }
    for(int i = 0; i < times; i++){
        temp = nextBinary(config, unionSnpCount);
    }
   // printf("num causal in this config is %d\n", temp);
    return config;
}

vector<int> PostCal::constructConfig(vector<int> input_causal_locs) {
   vector<int> config(totalSnpCount, 0);
   for ( int i = 0; i < input_causal_locs.size(); i++ ) {
      if ( input_causal_locs[i] >= 0 ) {
         config[input_causal_locs[i]] = 1;
      }
   }
   return config;
}


double PostCal::computeTotalLikelihoodGivenConfigs(vector<double>* stat, double sigma_g_squared) {
	printf("Input configs given\n");
    double sumLikelihood = 0;
    long int total_iteration = 0 ;
    int mycount = 0;
    printf("num total configs = %d\n", mycount);

    int unionSnpCount = all_snp_pos.size();
    cout << "Max Causal = " << maxCausalSNP << endl;
    cout << "Union Snp Count = " << unionSnpCount << endl;

    //clock_t start = clock();
    vector<int> configure;

    int chunksize;
    if(total_iteration < 1000){
	    if ( total_iteration < 10 ) {
               chunksize = total_iteration;
	    } else {
        chunksize = total_iteration/10;
	    }
    }
    else{
        chunksize = total_iteration/1000;
    }
    int curr_iter = 0;

    char *configs = NULL;
    size_t size_configs_file = 0;
    int status = safe_mmap_read_only(configsFile, &configs, &size_configs_file);
    if ( status < 0 ) {
      printf("mmap did not succeed\n");
      exit(1);
    }
    if ((num_configs * num_groups * sizeof(int16_t)) != size_configs_file) {
      printf("config file is not the expected size\n");
      exit(1);
    }

    //TODO fix this pragma
    //#pragma omp parallel for schedule(static,chunksize) private(configure)
    #pragma omp parallel for
    for(long int cidx = 0; cidx < num_configs; cidx++) {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	int16_t * input_causal_locs = (int16_t *) (configs + (cidx * num_groups * sizeof(int16_t))); //offset for (cidx)th config 

        vector<int> config(totalSnpCount, 0);
	int numCausal = 0;
        for ( int i = 0; i < num_groups; i++ ) {
          if ( input_causal_locs[i] >= 0 ) {
            config[input_causal_locs[i]] = 1;
	    numCausal += 1;
          }
        }
	//printVec(config);


        if ( numCausal == 0 ) { //if no causal, just update sum likelihood, nothing else should change
                vector<int> tempConfigure(totalSnpCount, 0);
                double tmp_likelihood = 0;

                if(haslowrank==true){
		mycount += 1;
		  mat tmpResultMatrixNN = statMatrixtTran * statMatrix;
                  double res = tmpResultMatrixNN(0,0);
                  double matDet = 1;
                  double lrl = (-res/2-sqrt(abs(matDet))); //lrl = low rank likelihood
                  tmp_likelihood = lrl + unionSnpCount * log(1-gamma);
                }
                else{
		  mycount += 1;
		  mat tmpResultMatrixNM = statMatrixtTran * invSigmaMatrix;
                  mat tmpResultMatrixNN = tmpResultMatrixNM * statMatrix;

                  double res = tmpResultMatrixNN(0,0);
                  double matDet = sigmaDet;
                  double ll = (-res/2-sqrt(abs(matDet)));
                  tmp_likelihood = ll + unionSnpCount * log(1-gamma);
                }  
		for ( int ss = 0; ss < num_of_studies; ss++ ) {
                  noCausal[ss] = addlogSpace(noCausal[ss], tmp_likelihood);
		}

                #pragma omp critical
                sumLikelihood = addlogSpace(sumLikelihood, tmp_likelihood);
		continue;
	}



	vector<int> causal_locs;
	//printf("causal locs are: ");
	int curr_study = 0;
	int curr_cum_num_snps = num_snps_all[curr_study];
	for ( int i = 0; i < num_groups; i++ ) {
	    int causal_snp_global_idx = input_causal_locs[i];
            if ( causal_snp_global_idx >= 0 ) {
	      //find which study this is associated with and find the snp's global position i.e. union pos
	      while (causal_snp_global_idx >= curr_cum_num_snps) {
                curr_study += 1;
		curr_cum_num_snps += num_snps_all[curr_study];
	      }
	      int union_pos = idx_to_union_pos_map[curr_study][causal_snp_global_idx - (curr_cum_num_snps - num_snps_all[curr_study])];
	      causal_locs.push_back(union_pos);
	    }    
	}
	//at this point, causal locs could have duplicates and is not sorted
	std::sort( causal_locs.begin(), causal_locs.end() );
        causal_locs.erase( unique( causal_locs.begin(), causal_locs.end() ), causal_locs.end() );
	//reset num causal to causal_locs.size(). num causal is now the unique number of causal snps in this config
	numCausal = (int)causal_locs.size();

	//printf("\n");


        int** causal_bool_per_study_for_config = new int*[num_of_studies]; //this implicitly mallocs so it must be freed
        int** causal_idx_per_study = new int*[num_of_studies]; //this implicitly mallocs so it must be freed
	for ( int i = 0; i < num_of_studies; i++ ) {
          causal_bool_per_study_for_config[i] = new int[numCausal];
          causal_idx_per_study[i] = new int[numCausal];
	}
	for ( int i = 0; i < num_of_studies; i++ ) {
	  memset(causal_bool_per_study_for_config[i], 0, numCausal * sizeof(int));
	  memset(causal_idx_per_study[i], 0, numCausal * sizeof(int));
	}

	/*printf("Start causal bool per study for config\n");
        for ( int i = 0; i < num_of_studies; i++ ) {
          for ( int j = 0; j < numCausal; j++ ) {
            printf("%d ", causal_bool_per_study_for_config[i][j]);
	  }
	  printf("\n");
	}
	printf("End causal bool per study for config\n");
	printf("Start causal idx per study for config\n");
        for ( int i = 0; i < num_of_studies; i++ ) {
          for ( int j = 0; j < numCausal; j++ ) {
            printf("%d ", causal_idx_per_study[i][j]);
	  }
	  printf("\n");
	}
	printf("End idx per study for config\n");
	*/

	int aux_idx = 0;
	while ( aux_idx < num_groups) { //find first not -1 entry
          if ( input_causal_locs[aux_idx] < 0 ) {
            aux_idx += 1;
	  } else {
            break;
	  }
	}
	curr_study = 0;
	curr_cum_num_snps = 0;
	for ( int i = 0; i < num_of_studies; i++ ) {
	    curr_cum_num_snps += num_snps_all[i];
            for ( int j = 0; j < numCausal; j++ ) {
		int loc_in_studyi = idx_to_snp_map[i][causal_locs[j]];
		//printf("loc in study %d is %d\n", i, loc_in_studyi);
		int studyi_offset = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
                if (loc_in_studyi >= 0) { //means it exists in study i
	          int global_idx = studyi_offset + loc_in_studyi;
		  if ( global_idx >= curr_cum_num_snps ) {
                    curr_study += 1; //none causal in ith study, skip to next study
		    //curr_cum_num_snps += num_snps_all[curr_study];
		    break;
		  }
		  else if ( input_causal_locs[aux_idx] == global_idx ) { //we have caught up to aux_idx location
		    aux_idx += 1; //need to increment aux idx to next loc
		    causal_bool_per_study_for_config[i][j] = 1;
		    causal_idx_per_study[i][j] = loc_in_studyi;
		    while ( aux_idx < num_groups) { //reset aux idx to next causal loc
                      if ( input_causal_locs[aux_idx] < 0 ) {
                        aux_idx += 1;
		      } else {
                        break;
		      }
		    }
		  }
		}
	    } 
	    if ( aux_idx == num_groups ) { //if we have incremented aux idx past the end, then we are all done
              break;
	    }
	}
	if ( aux_idx != num_groups) {
          printf("This did not work as expected\n");
	  exit(1);
	}
	
	/*
	printf("Start causal bool per study for config\n");
        for ( int i = 0; i < num_of_studies; i++ ) {
          for ( int j = 0; j < numCausal; j++ ) {
            printf("%d ", causal_bool_per_study_for_config[i][j]);
	  }
	  printf("\n");
	}
	printf("End causal bool per study for config\n");
	printf("Start causal idx per study for config\n");
        for ( int i = 0; i < num_of_studies; i++ ) {
          for ( int j = 0; j < numCausal; j++ ) {
            printf("%d ", causal_idx_per_study[i][j]);
	  }
	  printf("\n");
	}
	printf("End idx per study for config\n");
	*/

        double tmp_likelihood = 0;
	double just_ll = 0;
	double just_prior = 0;
        mat sigmaC = construct_diagC(config, numCausal, causal_idx_per_study, causal_bool_per_study_for_config);
        //printf("sigma C\n");
	//sigmaC.print(std::cout);

        if ( haslowrank == true ) {
           mycount += 1;
	   just_prior = 
           just_ll = lowrank_likelihood(config, stat, sigma_g_squared, sigmaC);
           just_prior =  log_prior(config, numCausal, causal_bool_per_study_for_config); 
           tmp_likelihood = just_ll + just_prior;
           }
        else {
           mycount += 1;
           just_ll = likelihood(config, stat, sigma_g_squared, sigmaC);
           just_prior = log_prior(config, numCausal, causal_bool_per_study_for_config);
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
                postValues[f] = addlogSpace(postValues[f], tmp_likelihood * config[f]);
		//if ( nextConfigure[f] == 1 ) {
                  //printf("updating index %d\n", f);
		  //printf("added in log space %f\n", tmp_likelihood);
		  //printf("post value for %d is %f\n", f, postValues[f]);
		//}
                //}
         }        
	 //printf("post values\n");
	 //for ( int f = 0; f < totalSnpCount; f++ ) {
         //  printf("%f ", postValues[f]);
	 //}
	 //printf("\n");

	 for(int i = 0; i < num_of_studies; ++i) {
           delete[] causal_bool_per_study_for_config[i];
           delete[] causal_idx_per_study[i];
         }
         //Free the array of pointers
         delete[] causal_bool_per_study_for_config;
         delete[] causal_idx_per_study;

        #pragma omp critical
        if(cidx % 1000 == 0 and cidx > curr_iter){
            cerr << "\r                                                                 \r" << (double) (cidx) / (double) num_configs * 100.0 << "%";
            curr_iter = cidx;
        }
	if(cidx % 1000 == 0 ){
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        //std::cout << "Time to eval config = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
	}
    }

    //cout << "\ncomputing likelihood of all configurations took  " << (float)(clock()-start)/CLOCKS_PER_SEC << "seconds.\n";

    printf("num total configs = %d\n", mycount);

    status = munmap(configs, size_configs_file);

    return(sumLikelihood);
}

double PostCal::computeTotalLikelihood(vector<double>* stat, double sigma_g_squared) {
    double sumLikelihood = 0;
    long int total_iteration = 0 ;
    int mycount = 0;
    printf("num total configs = %d\n", mycount);

    int unionSnpCount = all_snp_pos.size();
    for(long int i = 0; i <= maxCausalSNP; i++)
        //total_iteration = total_iteration + nCr(snpCount, i);
        total_iteration = total_iteration + nCr(unionSnpCount, i);
    cout << "Max Causal = " << maxCausalSNP << endl;
    cout << "Union Snp Count = " << unionSnpCount << endl;

    //clock_t start = clock();
    vector<int> configure;
    int num;

    int chunksize;
    if(total_iteration < 1000){
	    if ( total_iteration < 10 ) {
               chunksize = total_iteration;
	    } else {
        chunksize = total_iteration/10;
	    }
    }
    else{
        chunksize = total_iteration/1000;
    }
    int curr_iter = 0;

    
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


    #pragma omp parallel for schedule(static,chunksize) private(configure,num)
    for(long int iter = 0; iter < total_iteration; iter++) {
        if(iter%chunksize == 0){
            configure = findConfig(iter);
        }
        else{
            //num = nextBinary(configure, snpCount);
            num = nextBinary(configure, unionSnpCount);
        }


	int numCausal = std::accumulate(configure.begin(), configure.end(), 0);
	//if ((numCausal != 0) and (numCausal != maxCausalSNP)) {
	//	printf("%d is numCausal: ", numCausal);
	//	printf("%d is maxCausal: ", maxCausalSNP);
	//printf("Skipping initial configure: ");
	//printf("next configure to expand\n");
	//printVec(configure);
         // continue;
	//}
	//printf("Initial configure: ");
	//printVec(configure);
	//printf("numCausal = %d\n", numCausal);

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
		for ( int ss = 0; ss < num_of_studies; ss++ ) {
                  noCausal[ss] = addlogSpace(noCausal[ss], tmp_likelihood);
		}

                #pragma omp critical
                sumLikelihood = addlogSpace(sumLikelihood, tmp_likelihood);
		continue;
	}


	vector<int> causal_locs;
	//printf("causal locs are: ");
	for ( int i = 0; i < unionSnpCount; i++ ) {
            if ( configure[i] == 1 ) {
	        causal_locs.push_back(i);
		//printf("%d ", i);
	    }    
	}
	//printf("\n");
        
	
	//printf("total snp count is: %d\n", totalSnpCount);
        vector<int> startConfigure(totalSnpCount, 0);

	/*
	int** causal_idx_per_study = new int*[num_of_studies];
	int** causal_bool_per_study = new int*[num_of_studies];
	int** causal_bool_per_study_for_config = new int*[num_of_studies];
	for ( int i = 0; i < num_of_studies; i++ ) { //3 is hardcoded as max causal for now
	  causal_idx_per_study[i] = new int[3];
	  causal_bool_per_study[i] = new int[3];
	  causal_bool_per_study_for_config[i] = new int[3];
	}
	*/
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
		/*
          //int **causal_bool_per_study_for_config = new int*[num_of_studies];
	  int** causal_bool_per_study_for_config = (int**)malloc(num_of_studies * sizeof(int*));
	  if (causal_bool_per_study_for_config == nullptr) {printf("alloc error\n");}
	  for ( int j = 0; j < num_of_studies; j++ ) {
            //causal_bool_per_study_for_config[j] = new int[3];
	    causal_bool_per_study_for_config[j] = (int*)malloc(3 * sizeof(int));
	    if (causal_bool_per_study_for_config[j] == nullptr) {printf("alloc error\n");}
	  }
	  for ( int j = 0; j < num_of_studies; j++ ) {
	    memset(causal_bool_per_study_for_config[j], 0, 3 * sizeof(int));
	  }
	  */

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
	 /*
	 for ( int j = 0; j < num_of_studies; j++ ) {
           //delete[] causal_bool_per_study_for_config[j];
           free(causal_bool_per_study_for_config[j]);
	 }
	 //delete[] causal_bool_per_study_for_config;
	 free(causal_bool_per_study_for_config);
	 */
       }

	/*
	for(int i = 0; i < num_of_studies; ++i) {
           delete[] causal_idx_per_study[i];
           delete[] causal_bool_per_study[i];
           delete[] causal_bool_per_study_for_config[i];
         }
         //Free the array of pointers
         delete[] causal_idx_per_study;
         delete[] causal_bool_per_study;
         delete[] causal_bool_per_study_for_config;
	*/

        #pragma omp critical
        if(iter % 1000 == 0 and iter > curr_iter){
            cerr << "\r                                                                 \r" << (double) (iter) / (double) total_iteration * 100.0 << "%";
            curr_iter = iter;
        }
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
