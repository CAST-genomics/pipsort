# PIPSORT

## About

PIPSORT is a multi-ancestry fine-mapping tool. It takes in summary statistics and LD matrices from two studies, typically of two different ancestries, and outputs ancestry-specific PIPs as well as other probabilities of interest (see our manuscript linked below). PIPSORT was designed to distinguish GWAS signals that are shared across ancestries versus ancestry-specific. The variants do not have to be matched across ancestries, meaning it can handle variants that were tested in one but not both studies.

PIPSORT extends https://github.com/nlapier2/MsCAVIAR/

The majority of PIPSORT development was done on a fork. We have moved the code to this new repo. For prior commit history, please see https://github.com/TaraMirmira/MsCAVIAR

## Installation

Prerequisites:
- A recent version of C/C++ compiler supporting `C++11` standard
- `CMake` version `3.16` or above

Required libraries:
- GNU scientific library (GSL)
- BLAS and LAPACK
- C++ compiler

To compile:

```
git clone https://github.com/CAST-genomics/pipsort.git
cd pipsort/
mkdir build
cd build
cmake ..
make
```

This will generate the `PIPSORT` binary in the `build/` directory.

To install:

```
cmake --install . --prefix PREFIX
```

where `PREFIX` is a place you have write permissions. In most cases this will be your home directory, e.g. `$HOME`. If you install locally, make sure `$PREFIX/bin` is on your PATH.

Installation notes: in our development, we used GCC version 10.2.0, GSL version 2.5, and OpenBLAS version 0.3.27. Compiling PIPSORT generally takes less than a minute.

To test the install was successful, type `PIPSORT -h` which should show a help message.

## Quickstart

An small example command for running PIPSORT using files in this repository:

```bash
cd tests/small_example
PIPSORT -l ldfiles.txt \
	-z zfiles.txt \
	-m eur_afr_small_test_snp_map \
	-n 7000,7000 \
	-o test
```

This small example should take seconds to run. The input and output file formats, along with detailed usage, are described below.

## Usage

To run PIPSORT use the following command:

```
PIPSORT [OPTIONS] \
	-l <LDFILE> \
	-z <ZFile> \
	-m <snpMapFile> \
	-n <int,int> \
	-o <string>
```

### Required command line arguments

* **`-l <LDFile>`**: File containing paths to ld files (see input format below)
* **`-z <ZFile>`**: File containing paths to Z score files (see input format below) 
* **`-m <snpMapFile>`**: File mapping indexes to SNPs in each study (see input format below)
* **`-n <int,int>`**: Sample sizes (integers) of each study. e.g. 50,100. These should be the same order as was given by the `-l` and `-z` arguments
* **`-o <OUTPREFIX`**: Output prefix for PIPSORT output files (see outputs and formats below)

### Additional optional parameters

The most important optional parameters are:

* **`-c <int>`**: controls the maximum number of causal SNPs allowed at a locus; the default is 3. 
* **`-p <float>`**: sharing parameter giving the prior probability that a causal variant is shared across studies (default 0.75). See "Recommendations on parameter selection" below.

The parameters below can also be modified using command line options but we recommend sticking with the defaults:

* **`-t <float>`**: Sets the heterogeneity (t^2) across studies, default is 0.52
* **`-g <float>`**: Sets the prior of a SNP being causal (default 0.01)
* **`-s <float>`**: Sets the NCP variance for the smallest study, default is 5.2.

### Stochastic shotgun search

To improve the runtime burden, PIPSORT can be run with stochastic shotgun search to shrink the otherwise exhaustive search space. The command line option for this is:

* **`-q 1`** to use stochastic shotgun search, 0 (default) for exhaustive search 
* **`-x <int>`** to set the random seed. Default: 12345.

## Input file formats

### LD File and Z files

The `-l` and `-z` options (examples in `tests/example/`) take in text files that list the paths to the LD matrix and Z scores for the variants in the locus, respectively. The file paths should be given either as absolute file paths or relative paths to the current directory. Both files should contain one line per study (so two total lines).

The LD file is a tab separated file giving the pairwise LD (r, not r-squared) for each pair of variants. The number of rows and columns should match the total number of variants. For each study, the number of lines in the LD file should match the number of lines in the Z file. The number of lines across studies need not match to allow for the sets of variants to differ across studies. 

The Z files should contain two tab-separated columns. The first column should be the variant name and the second should be the Z-score for that variant. Variants must be in the same order as in the LD file for each study. See the files in the tests/example/ folder for examples of these file types.

### SNP map file

To accommodate different sets of variants across studies, PIPSORT needs a variant map (`-m` option). For two studies, this is a 3-column comma-separated file. The first column is the variant name. The second column is the index of the variant in the first study (or -1 if it is not present in that study). The third column is the index of the variant in the second study (or -1 if it is not present). We provide helper scripts for formatting summary and statistics and constructing this map in `utils/`. Additionally we provide a script (`get_pipsort_inputs.sh`) that takes in two GWAS summary statistic files in PLINK format and outputs the variant map and summary statistics formatted appropriately for use by PIPSORT.

## Recommendations on parameter selection

* Maximum number of causal variants (`-c`): The default is 3. PIPSORT runtime is highly dependent on this value. Running with a single causal variant (`-c 1`) will dramatically decrease run time but may decrease accuracy if there are multiple causal variants in a locus. For running with a large number of causal variants, we recommend using the `-q 1` option to use stochastic shotgun search. We do not recommend running with more than 3 total causal variants especially at loci with more than 1,000 input variants per study.

* Sharing parameter (`-p`): The default is 0.75, but can be modified by the user. The parameter is used in the prior to adjust the weight of configurations. Values closer to 1 will assign higher weight to configurations that model shared signals and values closer to 0 will assign higher weight to configurations that model study-specific signals. In practice, we recommend users try different values and compare across results. In our paper, we try both 0.25 and 0.75 on real data and compare PIP values across both sets of results.

## Output files

PIPSORT will output 6 files:

- `${OUTPREFIX}_study0_pips.txt`: ancestry-specific PIPs for the first study
- `${OUTPREFIX}_study0_set.txt`: all variants with ancestry-specific PIP >= 0.5
- `${OUTPREFIX}_study1_pips.txt`: ancestry-specific PIPs for the second study
- `${OUTPREFIX}_study1_set.txt`: all variants with ancestry-specific PIP >= 0.5
- `${OUTPREFIX}_nocausal.txt`: a single-column two line file with P(no causal in study 0) as the first number and P(no causal in study 1) as the second
- `${OUTPREFIX}_shared_pips.txt`: shared PIPs

## Larger example

An example for running PIPSORT can be found in `tests/example`. All necessary files are provided as well as expected output files. The example can be run with `bash run_example.sh`.

There are a few probabilities computed after PIPSORT runs: global PIPs and not shared PIPs. Scripts for computing these are provided in `utils/` and example usage is provided in `run_example.sh`.

With a single processor, this example takes about 15 minutes to run. With 64 processors, it takes less than 2 minutes.

## Example analysis workflows

A typical workflow for running PIPSORT consists of:

1. Obtaining GWAS summary statistics for two studies of interest (typically the same trait analyzed in two different ancestry groups). These may be either computed directly by the user if individual-level data is directly available or obtained from published summary statistics.

2. Identifying fine-mapping trait regions. Typically these consist of approximately 1 Mb windows centered at lead variants. See our manuscript below for how we defined trait-regions for fine-mapping with PIPSORT.

3. Extract Z scores and LD matrices for each study in each trait region. Z-scores can be computed the effect size divided by the standard error of each variant. LD matrices (r, not r-squared) can be either computed from available individual level data (preferred) or obtained from a published reference panel (e.g. 1000 Genomes) matched on ancestry to each study. Prior to this step, we recommend filtering variants that are rare (e.g. MAF<0.01) in both studies. **Note, to improve computational efficiency**, you can also optionally loosely filter variants with very insignificant p-values (e.g. P>0.1). While this may slightly decrease accuracy it is in some cases required to make running PIPSORT tractable. We typically apply PIPSORT to fine-map regions with up to 1,000 variants per study.

4. Run PIPSORT using the instructions provided above.

Examples of PIPSORT worfklows can be found in the `scripts/` directory here: https://github.com/TaraMirmira/pipsort_workflows.

## Important utility scripts

After running `PIPSORT` and obtaining output files `study0_post.txt`, `study1_post.txt`, `shared_pips.txt`, global PIPs and not-shared PIPs can be computed with the provided utility scripts like this:

```
python get_global_pips.py study0_post.txt study1_post.txt shared_pips.txt global_pips.txt
python get_not_shared_pips.py shared_pips.txt global_pips.txt not_shared_pips.txt
```

Additional key utility scripts:
- `extract_all_snps_rsid.py` can be used to construct the `snp_map` (see pipeline in example)

## Interpretation of probabilities output by PIPSORT

PIPSORT computes multiple variant-level probabilities (PIPs) in addition to a locus-level probability which can be used to identify signals that are shared vs. ancestry-specific. 

* Ancestry-specific PIPs (variant level): These are similar to PIPs output by traditional fine-mappers but provide a separate PIP for each of the two studies. Similar to other types of PIPs, they are noisy when power to detect signals is low (small sample size, small effect size). The difference in ancestry-specific PIPs can be used to help classify if a variant is likely shared vs. shows an ancestry-specific effect (AS-E). Based on our simulation studies, we have found a difference in PIPs of >=0.1 to be indicative of a potential AS-E. However, we note that sample size imbalance can also lead to large differences in ancestry-specific PIPs.

* Global PIP (variant level): This PIP is comparable to the PIP output by fine-mappers that do not compute ancestry-specific PIPs. A higher global PIP suggests higher confidence that a variant is a causal variant, and additional values output by PIPSORT can be used to further classify if the variant is likely shared vs. an AS-V or AS-E. We further note that for any particular variant, Global PIP = Ancestry 1 PIP + Ancestry 2 PIP - Shared PIP.

* Shared PIP (variant level): This value can be used to quantify evidence that a causal variant is shared. For shared causal variants that we are well-powered to detect, the shared PIP is often high even with sample size imbalances and low effect sizes. 

* Not-shared PIP (variant level): This value can be used to quantify evidence that a variant is causal but not shared. It is similar to the difference in ancestry-specific PIPs. Note: for any particular variant, Not-shared PIP = Global PIP-Shared PIP.

* P(no causal) (locus-level): Unlike the PIPs described above, which are variant-level, P(no causal) is a locus-level probability. This probability is most helpful for detecting when there are no signals at a locus for an ancestry. Low ancestry-specific PIPs and high P(no causal) are indicative of no causal signals. However, this outcome could also be a result of low power e.g. due to sample size, which should be taken into account. 

## Citation

If you use our tool, please cite our paper!
https://www.medrxiv.org/content/10.1101/2025.11.13.25339614v1

