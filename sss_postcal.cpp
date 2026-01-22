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
#include <cassert>

using namespace arma;

vector<vector<int>> PostCal::get_nbdplus(vector<int> curr_causal_locs) {
    //each vector will be the list of indices corresponding to causal snps
    vector<vector<int>> plusone_configs;
    int curr_num_causal = curr_causal_locs.size();

    if ( curr_num_causal >= maxCausalSNP) {
       return plusone_configs; //going to cap at max causal
    }
    
    vector<int> curr_config(unionSnpCount, 0);
    for ( int i = 0; i < curr_num_causal; i++ ) {
        curr_config[curr_causal_locs[i]] = 1;
    }

    int num_plus_configs = curr_config.size() - curr_num_causal;
    for ( int i = 0; i < curr_config.size(); i++ ) {
      if ( curr_config[i]  == 0 ) {
          vector<int> new_config;
	  new_config.push_back(i);
	  new_config.insert(new_config.end(), curr_causal_locs.begin(), curr_causal_locs.end());
	  std::sort(new_config.begin(), new_config.end()); //sort ascending
	  plusone_configs.push_back(new_config);
      }
    }
    if ( plusone_configs.size() != num_plus_configs ) {
       printf("plusone configs WRONG size\n");
    }
    return plusone_configs;
}

vector<vector<int>> PostCal::get_nbdminus(vector<int> curr_causal_locs) {
    //each vector will be the list of indices corresponding to causal snps
    vector<vector<int>> minusone_configs;
    int curr_num_causal = curr_causal_locs.size();

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


vector<vector<int>> PostCal::get_nbdzero(vector<int> curr_causal_locs) {
    //each vector will be the list of indices corresponding to causal snps
    vector<vector<int>> nbdzero_configs;
    int curr_num_causal = curr_causal_locs.size();

    vector<int> curr_config(unionSnpCount, 0);
    for ( int i = 0; i < curr_num_causal; i++ ) {
        curr_config[curr_causal_locs[i]] = 1;
    }

    int num_nbdzero_configs = curr_num_causal * (curr_config.size() - curr_num_causal);
    vector<vector<int>> minusone_configs = get_nbdminus(curr_causal_locs); 
    for ( int i = 0; i < curr_config.size(); i++ ) {
	if ( curr_config[i] == 0 ) {
	    for ( vector<int> v : minusone_configs ) {
		vector<int> new_config;
		new_config.push_back(i);
		new_config.insert(new_config.end(), v.begin(), v.end());
		std::sort(new_config.begin(), new_config.end()); //sort ascending
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
    //vector<int> configure(unionSnpCount, 0);
    vector<int> causal_locs;
    //causal_locs.push_back(330);
    //causal_locs.push_back(181);



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
    //nP = 4;
    //nP = 1;
    vector<int> thread_info(nP, 0);
    vector<int> thread_iters(nP, 0);
    vector<std::chrono::microseconds> thread_time(nP, std::chrono::microseconds(0));
    
    omp_set_num_threads(nP);
    printf("num threads is %d\n", nP); 
    int num_expansions = 0;

    std::chrono::microseconds total_time{0};

    double old_sum_lkl = 0;

    int total_iteration = 1000; //TODO for now
    for(int iter = 0; iter < total_iteration; iter++) {
        
	if ( causal_locs.size() == 0 ) {
           printf("curr causal = all zeros vector\n");
        } else {
           printf("curr causal: \n");
           printVec(causal_locs);
        }


	vector<vector<int>> nbdzero = get_nbdzero(causal_locs);
	vector<vector<int>> nbdminus = get_nbdminus(causal_locs);
	vector<vector<int>> nbdplus = get_nbdplus(causal_locs);
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
	//if (explored_set.find(causal_locs) != explored_set.end()) {
	//    make_updates = false;
        //}
        auto it = config_hashmap.find(causal_locs);
        if (it != config_hashmap.end()) {
            make_updates = false; 
        }
	
	int num_expansions_for_curr = 0;
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	double tmp_likelihood = expand_and_compute_lkl(causal_locs, make_updates, stat, sigma_g_squared, &num_expansions_for_curr);
	printf("lkl of curr = %lf\n", tmp_likelihood);
        std::chrono::steady_clock::time_point temp_end = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(temp_end - start).count();
        std::cout << "time for expand= " << ms << " ms\n";

	//if ( num_expansions_for_curr == 0 ) {
        //   assert(causal_locs.size() == 0);
	//}
	__atomic_add_fetch(&num_expansions, num_expansions_for_curr, 0);

	//if ( make_updates ) {
        //   explored_set.insert(causal_locs);
	//}
	
        int num_new_configs = 0;
	
#define DEBUG
	vector<double> loglkls(nbd.size(), 0);
        vector<int> not_done;
	printf("nbd size for iter %d is %ld\n", iter, nbd.size());
	//TODO chunksize, static, schedule (stuff to test)
        #pragma omp parallel for schedule(static, 1)
	for ( int i = 0; i < nbd.size(); i++ ) {
#ifdef DEBUG
	   int tnum = omp_get_thread_num();
	   thread_iters[tnum] += 1;
           std::chrono::steady_clock::time_point thread_start = std::chrono::steady_clock::now();
#endif
	   vector<int> v = nbd[i];
           auto it = config_hashmap.find(v);
           if (it != config_hashmap.end()) {
               loglkls[i] = it->second;  
           } else {
	       int l_num_expansions = 0;
	       double v_lkl = expand_and_compute_lkl(v, true, stat, sigma_g_squared, &l_num_expansions); //compute lkl but no update
	       loglkls[i] = v_lkl;
               #pragma omp critical
               {
               not_done.push_back(i);
               }

         
#ifdef DEBUG
	       thread_info[tnum] += l_num_expansions;
	       //if ( l_num_expansions == 0 ) {
               //   assert(v.size() == 0);
	       //}
	       __atomic_add_fetch(&num_expansions, l_num_expansions, 0); //TODO
	       //printf("num expansions for iter=%d, i=%d is %d\n",iter, i, l_num_expansions);
               std::chrono::steady_clock::time_point thread_end = std::chrono::steady_clock::now();
               thread_time[tnum] += std::chrono::duration_cast<std::chrono::microseconds>(thread_end - thread_start);
#endif
           }
	}
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        total_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start);
 
        //check break condition
        if ( not_done.size() == 0 ) {
            printf("hit break condition\n");
            break; //we have seen no new configs, end
        }
        //check convergence condition
        if ( iter >= 100 ) {
            if ( (1 - exp(old_sum_lkl - sss_sum_lkl)) <= 0.001 ) {
               printf("hit convergence condition\n");
               break; //we have seen no new configs, end
            }
        }

#ifdef DEBUG
        for ( int i : not_done ) {
           assert(i <= nbd.size());
           assert(i >= 0);
        }
#endif

	//do the hashmap updates outside the pragma
        for ( int i : not_done ) {
           // Zero allocation, zero copying
           const auto& v = nbd[i];
           config_hashmap[v] = loglkls[i];
        }


	//sampling

	double weight_zero = 0.0; size_t zero_sample = nbd.size(); //nbd.size() is an invalid index
	double weight_minus = 0.0; size_t minus_sample = nbd.size();
	double weight_plus = 0.0; size_t plus_sample = nbd.size();
        
        printf("num zero, num minus, num_plus = %d,%d,%d\n", num_zero, num_minus, num_plus);

	//std::random_device rd; // random
	//std::mt19937 gen(rd()); // random
	if ( num_zero != 0 ) {
            vector<double> probs;
            auto max_it = std::max_element(loglkls.begin(), loglkls.begin()+num_zero);
            double max_log = *max_it;
            for ( int ii = 0; ii < num_zero; ii++ ) {
               probs.push_back(exp(loglkls[ii]-max_log));
            }
	    std::discrete_distribution<size_t> dist(probs.begin(), probs.end());
	    zero_sample = dist(gen);
	    weight_zero = std::accumulate(probs.begin(), probs.end(), 0.0);
	}
	if ( num_minus != 0 ) {
            vector<double> probs;
            auto max_it = std::max_element(loglkls.begin()+num_zero, loglkls.begin()+num_zero+num_minus);
            double max_log = *max_it;
            for ( int ii = num_zero; ii < num_zero+num_minus; ii++ ) {
               probs.push_back(exp(loglkls[ii]-max_log));
            }
	    std::discrete_distribution<size_t> dist(probs.begin(), probs.end());
	    minus_sample = dist(gen);
	    weight_minus = std::accumulate(probs.begin(), probs.end(), 0.0);
	}
	if ( num_plus != 0 ) {
            vector<double> probs;
            auto max_it = std::max_element(loglkls.begin()+num_zero+num_minus, loglkls.end());
            double max_log = *max_it;
            for ( int ii = num_zero+num_minus; ii < loglkls.size(); ii++ ) {
               probs.push_back(exp(loglkls[ii]-max_log));
            }
	    std::discrete_distribution<size_t> dist(probs.begin(), probs.end());
	    plus_sample = dist(gen);
	    weight_plus = std::accumulate(probs.begin(), probs.end(), 0.0);
	}
        printf("weight zero = %f\n", weight_zero);
        printf("weight minus = %f\n", weight_minus);
        printf("weight plus = %f\n", weight_plus);

	std::discrete_distribution<size_t> dist({weight_zero, weight_minus, weight_plus});
	size_t idx = dist(gen);
        printf("picking idx = %ld\n", idx);
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
	
 	auto max_iter = std::max_element(loglkls.begin(), loglkls.end());
        int max_index = std::distance(loglkls.begin(), max_iter);
        printVec(nbd[max_index]);
	printf("max of probs = %lf\n", *max_iter);
	printf("final idx = %ld\n", final_idx);
	causal_locs = nbd[final_idx];
        printf("old sum lkl = %f\n", old_sum_lkl);
        printf("sss sum lkl = %f\n", sss_sum_lkl);
        old_sum_lkl = sss_sum_lkl;

    }

    omp_set_num_threads(1);

    std::cout << "Total time = " << total_time.count() << " µs\n";
    
    //cout << "\ncomputing likelihood of all configurations took  " << (float)(clock()-start)/CLOCKS_PER_SEC << "seconds.\n";

    for(int i = 0; i <= maxCausalSNP; i++) { //TODO what is this for, do I need to change it
        histValues[i] = exp(histValues[i]-sumLikelihood);
    }
    
    printf("explored %d configs\n", num_expansions);
    for ( int i = 0; i < thread_info.size(); i++ ) {
        printf("thread %d did %d expansions\n", i, thread_info[i]);
        printf("thread %d did %d iters\n", i, thread_iters[i]);
        printf("thread %d time = %lld\n", i, (long long)thread_time[i].count());
    }


    return(sss_sum_lkl);
}

double PostCal::fake_expand(vector<int> causal_locs) {
///*
  double sum = 0;
  //int n = 10000000; int reps = 1;
  int n = 1000000; int reps = 20;
  //int n = 2048; int reps = 25000;
  vector<double> x(n, 0.2); 
  vector<int> y(n, 0.1);
  vector<int> z(n, 0);
  for (int i = 1; i < n; i++ ) {
     y[i] += y[i-1];
     x[i] += x[i-1];
  }
  for ( int j = 0; j < reps; j++ ) {
     for ( int i = 0; i < n; i++ ) {
        z[i] += x[i] * y[i];
        x[i] += 0.01;
        y[i] += 0.01;
     }
  }
  return z[0];
//*/
  
/*
  int n = 256;
  arma::mat A = arma::randu<arma::mat>(n, n);
  arma::mat B = arma::randu<arma::mat>(n, n);
double ret = 0;
for ( int i = 0; i < 40; i++ ) {
  arma::mat C = A * B;
  ret = C(0,0);
}
  return ret;
*/

/*
size_t n = 256;
   std::random_device rd;
    std::mt19937 gen(654);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<std::vector<double>> B(n, std::vector<double>(n));
    std::vector<std::vector<double>> C(n, std::vector<double>(n));

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            A[i][j] = dist(gen);
            B[i][j] = dist(gen);
        }
    }
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double sum = 0;
            for (size_t k = 0; k < n; k++ ) {
               sum += A[i][k]*B[k][j];
            }
            C[i][j] = sum;
        }
    }

    return C[0][0];
*/

}

double PostCal::expand_and_compute_lkl(vector<int> causal_locs, bool make_updates, vector<double> * stat, double sigma_g_squared, int * l_num_expansions) {

	vector<int> configure(unionSnpCount, 0);
	
	int num_expansions = 0;

	int numCausal = causal_locs.size();
        for ( int i = 0; i < causal_locs.size(); i++ ) {
            configure[causal_locs[i]] = 1;
        }

	//auto it = config_hashmap.find(causal_locs);
	//if (it != config_hashmap.end()) {
    	//    return it->second;  // access the double value associated with the key
	//}

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
		if ( make_updates ) { 
		   #pragma omp critical
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
	  num_expansions += 1;
	//VecHash hasher;
	//std::size_t h = hasher(causal_locs);

	//std::cout << h << "\n";
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
        
	 if ( make_updates ) {
            #pragma omp critical
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
           #pragma omp critical
	   if ( sharedCausal ) {
             sharedPips[causal_locs[g]] =  addlogSpace(sharedPips[causal_locs[g]], tmp_likelihood);
	     sharedLL[causal_locs[g]] = addlogSpace(sharedLL[causal_locs[g]], just_ll);
	   } else {
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
	//config_hashmap[causal_locs] = max_lkl;

	for (int i = 0; i < num_of_studies; i++) {
    		delete[] causal_idx_per_study[i];
    		delete[] causal_bool_per_study[i];
    		delete[] causal_bool_per_study_for_config[i];
	}
	delete[] causal_idx_per_study;
	delete[] causal_bool_per_study;
	delete[] causal_bool_per_study_for_config;

	*l_num_expansions = num_expansions;

	return max_lkl;
}

