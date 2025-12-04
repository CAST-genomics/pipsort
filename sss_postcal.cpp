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

vector<vector<int>> PostCal::get_nbdplus(vector<int> curr_config) {
    //each vector will be the list of indices corresponding to causal snps
    vector<vector<int>> plusone_configs;
    int curr_num_causal = 0;
    vector<int> curr_causal_locs;
    for ( int i = 0; i < curr_config.size(); i++ ) {
        if ( curr_config[i] == 1 ) {
           curr_num_causal += 1;
           curr_causal_locs.push_back(i);
	}
    }
    if ( curr_num_causal >= 10 ) {
       return plusone_configs; //going to cap at max 10 causal variants
    }
    int num_plus_configs = curr_config.size() - curr_num_causal;
    for ( int i = 0; i < curr_config.size(); i++ ) {
      if ( curr_config[i]  == 0 ) {
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
        if ( curr_config[i] == 1 ) {
           curr_num_causal += 1;
           curr_causal_locs.push_back(i);
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
        if ( curr_config[i] == 1 ) {
           curr_num_causal += 1;
           curr_causal_locs.push_back(i);
	}
    }
    int num_nbdzero_configs = curr_num_causal * (curr_config.size() - curr_num_causal);
    vector<vector<int>> minusone_configs = get_nbdminus(curr_config); 
    for ( int i = 0; i < curr_config.size(); i++ ) {
	if ( curr_config[i] == 0 ) {
	    for ( vector<int> v : minusone_configs ) {
		vector<int> new_config;
		new_config.push_back(i);
		new_config.insert(new_config.end(), v.begin(), v.end());
		nbdzero_configs.push_back(new_config);
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


    //clock_t start = clock();
    //TODO initialized as this for now
    vector<int> configure(unionSnpCount, 0);



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


    std::mt19937 gen(12345); //deterministic
    
    int nP = omp_get_num_procs();
    
    omp_set_num_threads(nP);
    printf("num threads is %d\n", nP);

    //TODO move the pragma somewhere else
    //#pragma omp parallel for schedule(static,chunksize) private(configure,num)
    int total_iteration = 3000; //TODO for now
    #pragma omp parallel for
    for(int iter = 0; iter < total_iteration; iter++) {

	int numCausal = 0;
	vector<int> causal_locs;
    	for ( int i = 0; i < configure.size(); i++ ) {
            if ( configure[i] == 1 ) {
               numCausal += 1;
               causal_locs.push_back(i);
            }
        }
	if ( numCausal == 0 ) {
            causal_locs.push_back(static_cast<int>(configure.size())); //repr for all zeros config
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
	//check if we have already explored this config
	if (explored_set.find(causal_locs) != explored_set.end()) {
	    make_updates = false;
        }
	
	double tmp_likelihood = expand_and_compute_lkl(configure, make_updates, stat, sigma_g_squared);

	if ( make_updates ) {
           explored_set.insert(causal_locs);
	}

	vector<double> probs(nbd.size(), 0);
	//TODO pragma here
	//change from for each to just for, change push_back to probs[i] = X
	for ( int i = 0; i < nbd.size(); i++ ) {
	   vector<int> v = nbd[i];
	   double v_lkl = expand_and_compute_lkl(v, false, stat, sigma_g_squared); //compute lkl but no update
	   probs[i] = v_lkl;
	}

	//sampling

	double weight_zero = 0.0; size_t zero_sample = nbd.size(); //nbd.size() is an invalid index
	double weight_minus = 0.0; size_t minus_sample = nbd.size();
	double weight_plus = 0.0; size_t plus_sample = nbd.size();

	//std::random_device rd; // random
	//std::mt19937 gen(rd()); // random
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
	if ( num_plus != 0 ) {
	    std::discrete_distribution<size_t> dist(probs.begin()+num_zero+num_minus, probs.end());
	    plus_sample = dist(gen);
	    weight_plus = std::accumulate(probs.begin()+num_zero+num_minus, probs.end(), 0.0);
	}

	std::discrete_distribution<size_t> dist({weight_zero, weight_minus, weight_plus});
	size_t idx = dist(gen);
	size_t final_idx;
	switch (idx) {
    	    case 0:
                final_idx = zero_sample;
        	break;
    	    case 1:
                final_idx = minus_sample + num_zero;
                break;
            case 2:
                final_idx = plus_sample + num_zero + num_minus;
                break;
        }
	vector<int> next_causal_locs = nbd[final_idx];
	std::fill(configure.begin(), configure.end(), 0);
	for ( int causal_idx : next_causal_locs ) {
            configure[causal_idx] = 1;
	}

    }

    omp_set_num_threads(1);

    
    //cout << "\ncomputing likelihood of all configurations took  " << (float)(clock()-start)/CLOCKS_PER_SEC << "seconds.\n";

    for(int i = 0; i <= maxCausalSNP; i++) { //TODO what is this for, do I need to change it
        histValues[i] = exp(histValues[i]-sumLikelihood);
    }


    return(sss_sum_lkl);
}


double PostCal::expand_and_compute_lkl(vector<int> configure, bool make_updates, vector<double> * stat, double sigma_g_squared) {

	int numCausal = 0;
        vector<int> causal_locs;
        for ( int i = 0; i < configure.size(); i++ ) {
            if ( configure[i] == 1 ) {
               numCausal += 1;
               causal_locs.push_back(i);
            }
        }

	/*auto it = config_hashmap.find(causal_locs);
	if (it != config_hashmap.end()) {
    	    return it->second;  // access the double value associated with the key
	}*/

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

		//TODO where to put the pragma
                #pragma omp critical
		if ( make_updates ) { 
                   sss_sum_lkl = addlogSpace(sss_sum_lkl, tmp_likelihood);
		   //emplace won't make a copy of causal_locs
		   //config_hashmap.emplace(causal_locs, tmp_likelihood);
		}
		return tmp_likelihood;
	}


	
	//printf("total snp count is: %d\n", totalSnpCount);
        vector<int> startConfigure(totalSnpCount, 0);


	int** causal_idx_per_study = new int*[num_of_studies];
	int** causal_bool_per_study = new int*[num_of_studies];
	int** causal_bool_per_study_for_config = new int*[num_of_studies];
	for (int i = 0; i < num_of_studies; i++) {
    		causal_idx_per_study[i] = new int[numCausal]();   // zero-initialize
    		causal_bool_per_study[i] = new int[numCausal]();   // zero-initialize
    		causal_bool_per_study_for_config[i] = new int[numCausal]();   // zero-initialize
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

	double max_lkl = 0.0;

	for ( int i = 0; i < total_num_additional_configs; i++ ) {
	  //int pid = omp_get_thread_num();
	  //if ( ( pid < 0 ) || ( pid >= nP ) ) { printf("BAD PID\n"); }
	  for ( int j = 0; j < num_of_studies; j++ ) {
             // std::fill(causal_bool_per_study_for_config[i], causal_bool_per_study_for_config[i] + numCausal, 0); <-- might be safer?
	    memset(causal_bool_per_study_for_config[j], 0, numCausal * sizeof(int));
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
	  //printf("next config to eval\n");
	  //printVec(nextConfigure);
          double tmp_likelihood = 0;
	  double just_ll = 0;
	  double just_prior = 0;
          mat sigmaC = construct_diagC(nextConfigure, numCausal, causal_idx_per_study, causal_bool_per_study_for_config);

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
          
	  //update max_lkl
	  if (abs(tmp_likelihood) > abs(max_lkl) ) {
             max_lkl = tmp_likelihood;
	  }
        
         #pragma omp critical
	 if ( make_updates ) {
            sss_sum_lkl = addlogSpace(sss_sum_lkl, tmp_likelihood);
	 }

         if ( make_updates ) {
         for ( int w = 0; w < num_of_studies; w++ ) {
	   bool allZero = true;
           for ( int v = 0; v < numCausal; v++ ) {
             if (causal_bool_per_study_for_config[w][v] == 1) {
	       allZero = false;
	       break;
	     }
	   }
	   if ( allZero ) {
                #pragma omp critical
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
                #pragma omp critical
             sharedPips[causal_locs[g]] =  addlogSpace(sharedPips[causal_locs[g]], tmp_likelihood);
	     sharedLL[causal_locs[g]] = addlogSpace(sharedLL[causal_locs[g]], just_ll);
	   } else {
                #pragma omp critical
             notSharedLL[causal_locs[g]] = addlogSpace(notSharedLL[causal_locs[g]], just_ll);
	   }
	 }

         for(int f = 0; f < totalSnpCount; f++) {
            //for(int k = 0; k < num_of_studies; k++){
                #pragma omp critical
                postValues[f] = addlogSpace(postValues[f], tmp_likelihood * nextConfigure[f]);
         }        
	 }
       }
	//#pragma omp critical
	//config_hashmap[causal_locs] = max_lkl;

	for (int i = 0; i < num_of_studies; i++) {
    		delete[] causal_idx_per_study[i];
    		delete[] causal_bool_per_study[i];
    		delete[] causal_bool_per_study_for_config[i];
	}
	delete[] causal_idx_per_study;
	delete[] causal_bool_per_study;
	delete[] causal_bool_per_study_for_config;

	return max_lkl;
}

