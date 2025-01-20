import os
import sys
import pandas as pd
import numpy as np


n = len(sys.argv)
assert(n == 5)

eur_post_file = sys.argv[1]
afr_post_file = sys.argv[2]
shared_pips_file = sys.argv[3]
outfile = sys.argv[4]

shared_pips = pd.read_csv(shared_pips_file, sep="\t")
eur_post = pd.read_csv(eur_post_file, sep="\t")
afr_post = pd.read_csv(afr_post_file, sep="\t")
merged = eur_post.merge(afr_post, how="outer", on="SNP_ID")
assert(merged.shape[0] == shared_pips.shape[0])
merged = merged.merge(shared_pips, how='inner', on='SNP_ID')
merged["Prob_in_pCausalSet_x"] = merged["Prob_in_pCausalSet_x"].fillna(0)
merged["Prob_in_pCausalSet_y"] = merged["Prob_in_pCausalSet_y"].fillna(0)
merged["global_pips"] = merged["Prob_in_pCausalSet_x"] + merged["Prob_in_pCausalSet_y"] - merged['shared_pip']

def less_than_or_close_to_1(x):
    return x < 1.0 or np.isclose(x, 1.0)


#assert((merged['global_pips']<=1).all() or np.allclose(merged['global_pips'], 1))

assert(((merged['global_pips']<=1.0)|(np.isclose(merged['global_pips'],1.0,atol=1e-06))).all())
assert(((merged['global_pips']>=0)|(np.isclose(merged['global_pips'],0,atol=1e-06))).all())
merged.drop(["Prob_in_pCausalSet_x", "Prob_in_pCausalSet_y", 'shared_pip'], inplace=True, axis=1)
print(merged)
merged.to_csv(outfile, sep="\t", index=False)
