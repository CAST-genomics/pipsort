import pandas as pd


import sys
 
# total arguments
n = len(sys.argv)
assert(n == 3)
 
infile = sys.argv[1]
outfile = sys.argv[2]

eur_gwas = pd.read_csv(infile, sep="\t")

cols_to_drop = ["T_STAT","P","TEST","OBS_CT","A1"] #this is what should be there to drop
actual_cols_to_drop = [] #this is what is actually there
for c in cols_to_drop:
    if c in eur_gwas.columns:
        actual_cols_to_drop.append(c)


eur_gwas.drop(actual_cols_to_drop,axis=1,inplace=True)

eur_gwas.rename(columns={'#CHROM':'CHROMOSOME','POS':'BP', 'REF': 'REF_ALLELE', 'ALT': 'ALT_ALLELE', 'ID':'SNP_ID'}, inplace=True)

eur_gwas.to_csv(outfile, sep="\t", index=False, na_rep='NA')

eur_gwas = pd.read_csv(outfile, sep="\t")

eur_gwas.dropna(axis=0, inplace=True)

eur_gwas.to_csv(outfile, sep="\t", index=False, na_rep='NA')

