import pandas as pd
import numpy as np
import sys
import os

n = len(sys.argv)
assert(n == 4)

shared_pips_file = sys.argv[1]
global_pips_file = sys.argv[2]
outfile = sys.argv[3]
#eur_gwas_file = sys.argv[3]
#afr_gwas_file = sys.argv[4]

shared_pips = pd.read_csv(shared_pips_file, sep="\t")
global_pips = pd.read_csv(global_pips_file, sep="\t")
#eur_gwas = pd.read_csv(eur_gwas_file, sep="\t")
#afr_gwas = pd.read_csv(afr_gwas_file, sep="\t")

#all_snps = eur_gwas.merge(afr_gwas, on=["BP", "SNP_ID"], how="outer")
#all_snps = all_snps[["BP", "SNP_ID"]]
#print(all_snps)

#shared_pips = shared_pips.merge(all_snps, left_on="POS", right_on="BP", how="left")
#shared_pips.drop(["POS","BP"], inplace=True, axis=1)
#shared_pips.to_csv(shared_pips_file, header=True, index=False, sep="\t")
assert(shared_pips.shape[0] == global_pips.shape[0])
not_shared = global_pips.merge(shared_pips, on="SNP_ID", how="inner")
print(not_shared)
not_shared['not_shared_pip'] = not_shared['global_pips'] - not_shared['shared_pip']
assert((not_shared['not_shared_pip']<=1).all())
assert(( (not_shared['not_shared_pip']>=0) | (np.isclose(not_shared['not_shared_pip'],0,atol=1e-06)) ).all())
assert((not_shared['shared_pip']<=1).all())
assert(( (not_shared['shared_pip']>=0) | (np.isclose(not_shared['shared_pip'],0,atol=1e-06)) ).all())
not_shared.drop(['global_pips', 'shared_pip'], inplace=True, axis=1)
not_shared.to_csv(outfile, header=True, index=False, sep="\t")
