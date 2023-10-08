../../MsCAVIAR -l ldfiles.txt -z zfiles.txt -m eur_afr_small_test_snp_map -n 7000,7000 -b all_configs_int16 -d 72 -e 5 -o mscaviar_results


if cmp -s "mscaviar_results_study0_post.txt" "expected_study0_post.txt"; then
    printf 'study 0 post: SUCCESS\n'
else
    printf 'study 0 post: FAILURE\n'
    exit 1
fi

if cmp -s "mscaviar_results_study1_post.txt" "expected_study1_post.txt"; then
    printf 'study 1 post: SUCCESS\n'
else
    printf 'study 1 post: FAILURE\n'
    exit 1
fi

if cmp -s "mscaviar_results_study0_set.txt" "expected_study0_set.txt"; then
    printf 'study 0 set: SUCCESS\n'
else
    printf 'study 0 set: FAILURE\n'
    exit 1
fi

if cmp -s "mscaviar_results_study1_set.txt" "expected_study1_set.txt"; then
    printf 'study 1 set: SUCCESS\n'
else
    printf 'study 1 set: FAILURE\n'
    exit 1
fi

if cmp -s "mscaviar_results_nocausal.txt" "expected_nocausal.txt"; then
    printf 'no causal: SUCCESS\n'
else
    printf 'no causal: FAILURE\n'
    exit 1
fi

if cmp -s "mscaviar_results_shared_pips.txt" "expected_shared_pips.txt"; then
    printf 'shared pips: SUCCESS\n'
else
    printf 'shared pips: FAILURE\n'
    exit 1
fi

#rm mscaviar_results_study0_set.txt
#rm mscaviar_results_study1_set.txt
#rm mscaviar_results_study0_post.txt
#rm mscaviar_results_study1_post.txt
#rm mscaviar_results_shared_pips.txt
#rm mscaviar_results_nocausal.txt
