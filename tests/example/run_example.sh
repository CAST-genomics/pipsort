../../PIPSORT -c 2 -l ldfiles.txt -z zfiles.txt -m snp_map -n 334324,6771 -p 0.25 -o pipsort_results

utils=../../utils
python $utils/get_global_pips.py pipsort_results_study0_post.txt pipsort_results_study1_post.txt pipsort_results_shared_pips.txt global_pips.txt
python $utils/get_not_shared_pips.py pipsort_results_shared_pips.txt global_pips.txt not_shared_pips.txt
