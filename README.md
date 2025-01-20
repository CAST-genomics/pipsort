Extends https://github.com/nlapier2/MsCAVIAR/

### Installation


```
git clone https://github.com/CAST-genomics/pipsort.git
cd pipsort/
make
```


### Important command line options

The **-l** and **-z** arguments, as can be seen in the example/ folder, list the files with the LD matrix and Z scores for the SNPs in the locus, respectively. The file paths should be given either as absolute file paths or relative paths to the current directory. The LD files for a locus with M SNPs should contain M lines, each with M space-separated numbers. The _i_ th number in the _j_ th line should be the LD between SNPs _i_ and _j_. 

The Z files should contain M lines, each with two tab-separated columns. The first column should be the SNP name and the second should be the Z-score for that SNP. See the files in the example/ folder for examples of these file types. Please make sure the exact same set of SNPs are present for all files and in the same order. Also, make sure the studies are in the same order for -l and -z.

**-c** controls the maximum number of causal SNPs allowed at a locus; the default is 3. **-r** ("rho") determines the posterior probability threshold for MsCAVIAR; in other words, MsCAVIAR will return a set of SNPs that, with rho% probability, contains all causal SNPs. The default is 0.95. **-t** controls the heterogeneity between the studies; in other words, the variance in the effect sizes of causal SNPs not accounted for by sample size imbalance. The default is 0.52.

The **-n** argument specifies the population size for each study, comma-separated. Again, this should be the same number of studies and in the same order as was given by the -l and -z arguments. 

The **-o** argument is simply the output name prefix for the MsCAVIAR output files.

### Other command line options (from MsCAVIAR)

The other command line options are listed below. We do not recommend changing -g or -s for most users.


**-f** allows MsCAVIAR to output a _hist.txt_ file that includes the likelihood of the true set containing 0, 1, 2, ... up to the maximum number of causal SNPs that the user sets using the -c option. We recommend using this option when the user is unsure about the number of causal snps in the study, and then adjust -c according to the likelihoods in the histogram file. For example, if the histogram file contains 5 numbers: "1.21973e-150"  "2.54771e-09"  "0.614328"  "0.334465"  "0.051207", it means that the causal set most likely contains 2 snps, with L(c=2) = 0.61. But since L(c=3)=0.33 indicates that c=3 is also relatively likely, we recommend the user to rerun the program using -c 2 or -c 3 as the parameters to get more accurate results.

**-a** sets the threshold for cut-off when we select the final causal set: if a snp has a probability less than the threshold, then we do not include it in the set. Note that a high cutoff threshold may result in the set having a confidence less than 95%.

```
-r RHO, --rho-prob=RHO     set $rho$ probability (default 0.95)
-g GAMMA, --gamma      set $gamma$ the prior of a SNP being causal (default 0.01)
-c causal          set the maximum number of causal SNPs (default 3)
-f 1               to out the probaility of different number of causal SNP
-t TAU_SQR, --tau_sqr=TAU_SQR  set the heterogeneity (t^2) across studies, default is 0.52
-s SIGMA_G_SQR, --sigma_g_squared=SIGMA_G_SQR    set the NCP variance for the smallest study, default is 5.2
-a THRESHOLD    set the threshold for cut-off when we select the final causal set, default is 0

```

### Example

