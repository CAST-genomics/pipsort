# PIPSORT

### About

PIPSORT is a multi-ancestry fine-mapping tool. It takes in summary statistics and LD matrices from studies of two ancestries and outputs ancestry-specific PIPs as well as other probabilities of interest (see our manuscript linked below). PIPSORT was designed to distinguish GWAS signals that are shared across ancestries versus ancestry-specific. The variants do not have to be matched across ancestries. 

PIPSORT extends https://github.com/nlapier2/MsCAVIAR/

The majority of PIPSORT development was done on a fork. We have moved the code to this new repo. For prior commit history, please see https://github.com/TaraMirmira/MsCAVIAR

### Installation

Required libraries:
- GNU scientific library
- BLAS and LAPACK
- C++ compiler

```
git clone https://github.com/CAST-genomics/pipsort.git
cd pipsort/
make
```


### Required command line options

The **-l** and **-z** arguments, as can be seen in the tests/example/ folder, list the files with the LD matrix and Z scores for the SNPs in the locus, respectively. The file paths should be given either as absolute file paths or relative paths to the current directory.

For each study, the number of lines in the LD file should match the number of lines in the Z file. The number of lines across studies need not match to allow for the sets of variants to differ across studies. 

The Z files should contain two tab-separated columns. The first column should be the variant name and the second should be the Z-score for that variant. See the files in the tests/example/ folder for examples of these file types.

**-m** To accommodate different sets of variants across studies, PIPSORT needs a variant map. For two studies, this is a 3-column comma-separated file. The first column is the variant name. The second column is the index of the variant in the first study (or -1 if it is not present in that study). The third column is the index of the variant in the second study (or -1 if it is not present). We provide helper scripts for formatting summary and statistics and constructing this map in `utils/`. Additionally we provide a script (`get_pipsort_inputs.sh`) that takes in two GWAS summary statistic files in PLINK format and outputs the variant map and summary statistics formatted appropriately for use by PIPSORT.

The **-n** argument specifies the population size for each study, comma-separated. Again, this should be the same number of studies and in the same order as was given by the -l and -z arguments. 

The **-o** argument is the output name prefix for the PIPSORT output files.

### Important command line options

**-c** controls the maximum number of causal SNPs allowed at a locus; the default is 3. 

**-sp** is the sharing parameter. The default is 0.75, but can be modified by the user. The parameter is used in the prior to adjust the weight of configurations. Values closer to 1 will assign higher weight to configurations that model shared signals and values closer to 0 will assign higher weight to configurations that model study-specific signals. In practice, we recommend users try different values and compare across results. In our paper, we try both 0.25 and 0.75 on real data and compare PIP values across both sets of results.  

### Other command line options (from MsCAVIAR)

The other command line options are listed below. We do not recommend changing these for most users.

**-t** controls the heterogeneity between the studies; in other words, the variance in the effect sizes of causal SNPs not accounted for by sample size imbalance. The default is 0.52.

```
-g GAMMA, --gamma      set $gamma$ the prior of a SNP being causal (default 0.01)
-t TAU_SQR, --tau_sqr=TAU_SQR  set the heterogeneity (t^2) across studies, default is 0.52
-s SIGMA_G_SQR, --sigma_g_squared=SIGMA_G_SQR    set the NCP variance for the smallest study, default is 5.2

```

### Example

An example for running PIPSORT can be found in `tests/example`. All necessary files are provided as well as expected output files. The example can be run with `run_example.sh`.

There are a few probabilities computed after PIPSORT runs: global PIPs and not shared PIPs. Scripts for computing these are provided in `utils/` and example usage is provided in `run_example.sh`.

### Stochastic Shotgun Search

To improve the runtime burden, PIPSORT can be run with stochastic shotgun search to shrink the otherwise exhaustive search space. The command line option for this is:

**-q** 1 to use stochastic shotgun search, 0 (default) for exhaustive search 

For reproduciblity, we set a seed for stochastic shotgun search. We will change this to a user-specified parameter, but for now this can be removed/modified in line 138 of `sss_postcal.cpp`: https://github.com/CAST-genomics/pipsort/blob/main/sss_postcal.cpp#L138

### Citation

If you use our tool, please cite our paper!
https://www.medrxiv.org/content/10.1101/2025.11.13.25339614v1

