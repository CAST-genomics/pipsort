import pandas as pd
import numpy as np
import sys

gwas_file = sys.argv[1]
gwas = pd.read_csv(gwas_file, sep="\t")

#cols of importance are POS and ID

gwas.sort_values("POS", inplace=True, ignore_index=True)

orig_idx = np.arange(gwas.shape[0], dtype=int)
curr_pos = gwas["POS"][0]
curr_idx = 0
for i in range(1, gwas.shape[0]):
    next_pos = gwas["POS"][i]
    if (next_pos != curr_pos):
        curr_idx += 1
    orig_idx[i] = curr_idx
    curr_pos = next_pos
gwas["orig_idx"] = orig_idx
gwas.sort_values(["orig_idx","ID"], inplace=True, ignore_index=True)
gwas.drop(columns=["orig_idx"],inplace=True)
gwas.to_csv(gwas_file, sep="\t", index=False)
