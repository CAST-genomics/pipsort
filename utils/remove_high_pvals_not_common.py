import pandas as pd
import sys

n = len(sys.argv)

eur_infile = sys.argv[1]
afr_infile = sys.argv[2]
fsl1 = float(sys.argv[3])
fsl2 = fsl1
if n == 5:
    fsl2 = float(sys.argv[4])


df_eur = pd.read_csv(eur_infile, sep="\t")
df_afr = pd.read_csv(afr_infile, sep="\t")


df_eur_sig = df_eur.merge(df_afr, on="ID", how="left")
df_eur_sig = df_eur_sig[((df_eur_sig['P_x'] != pd.NA) & (df_eur_sig['P_x'] < fsl1)) | ((df_eur_sig['P_y'] != pd.NA) & (df_eur_sig['P_y'] < fsl2))]

df_afr_sig = df_afr.merge(df_eur, on="ID", how="left")
df_afr_sig = df_afr_sig[((df_afr_sig['P_x'] != pd.NA) & (df_afr_sig['P_x'] < fsl1)) | ((df_afr_sig['P_y'] != pd.NA) & (df_afr_sig['P_y'] < fsl1))]

df_eur_sig_final = df_eur.merge(df_eur_sig['ID'], on="ID", how="inner")
df_afr_sig_final = df_afr.merge(df_afr_sig['ID'], on="ID", how="inner")

df_eur_sig_final.to_csv("s1_gwas_temp.txt", sep="\t", index=False)
df_afr_sig_final.to_csv("s2_gwas_temp.txt", sep="\t", index=False)

df_eur_sig_final["ID"].to_csv("s1_snps.txt", header=False, index=False)
df_afr_sig_final["ID"].to_csv("s2_snps.txt", header=False, index=False)

