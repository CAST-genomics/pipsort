s1_gwas=$1
s2_gwas=$2


fsl_s1=0.0001 #filter significance level, default is 0.0001
fsl_s2=0.0001 #filter significance level, default is 0.0001


popsize1=$(awk 'NR==1 {for (i=1; i<=NF; i++) if ($i == "OBS_CT") col=i} NR==2 {print $col}' $s1_gwas)
echo $popsize1
popsize2=$(awk 'NR==1 {for (i=1; i<=NF; i++) if ($i == "OBS_CT") col=i} NR==2 {print $col}' $s2_gwas)
echo $popsize2


#needed to keep address snps with same pos but different ids
python sort_gwas_results.py $s1_gwas
python sort_gwas_results.py $s2_gwas


#exclude high p vals, outputs s1_snps.txt and s2_snps.txt
python remove_high_pvals_not_common.py $s1_gwas $s2_gwas $fsl_s1 $fsl_s2

python get_gwas_data_in_pipsort_style.py s1_gwas_temp.txt s1_gwas.txt
python get_gwas_data_in_pipsort_style.py s2_gwas_temp.txt s2_gwas.txt


python format_sumstats.py --infile s1_gwas.txt --outfile temp_f_s1.txt --chromosome CHROMOSOME --bp BP --snp_id SNP_ID --ref_allele REF_ALLELE --alt_allele ALT_ALLELE --beta BETA --se SE
python format_sumstats.py --infile s2_gwas.txt --outfile temp_f_s2.txt --chromosome CHROMOSOME --bp BP --snp_id SNP_ID --ref_allele REF_ALLELE --alt_allele ALT_ALLELE --beta BETA --se SE

#drop any duplicate rsid and zscore pairs, happens because of STR duplicates
awk '!seen[$0]++' temp_f_s1.txt > f_s1.txt
awk '!seen[$0]++' temp_f_s2.txt > f_s2.txt
rm temp_f_s1.txt
rm temp_f_s2.txt

python format_zscores.py f_s1.txt s1_final.zscores
python format_zscores.py f_s2.txt s2_final.zscores


echo $'f_s1.txt\nf_s2.txt' > processedfiles.txt
python extract_all_snps_rsid.py processedfiles.txt snp_map


rm f_s1.txt
rm f_s2.txt
rm processedfiles.txt
