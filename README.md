# HMPaint
**HMPaint** is an efficient tool for local ancestry inference (LAI) coded in C++. It improves **d-PBWT** algorithm to find K longest matches at each position, and uses the **Hash Map** strategy to implement the forward and backward algorithm in the Hidden Markov Model (HMM) because of the sparsity of haplotype matches. HMPaint incorporates the function for efficiently calculating [Linkage Disequilibrium of Ancestry (LDA), LDA score (LDAS)](https://github.com/YaolingYang/LDAandLDAscore) and [Ancestry Anomaly Score (AAS)](https://github.com/danjlawson/ms_paper) for understanding the population structure and evolution.

# Installation

The main code is in **devel/hashmap.cpp**.

You should load the [Armadillo](https://arma.sourceforge.net/download.html) library, and also have ["gzstream.h" and "gzstream.C"](https://www.cs.unc.edu/Research/compgeom/gzstream/) in your directory. 

To prepare the input files (phase format) for **HMPaint**, you should also get [PBWT](https://github.com/richarddurbin/pbwt) installed, which converts Variant Call Format (VCF) to phase format by the following command:

``
pbwt -readVcfGT XXX.vcf -writePhase XXX.phase
``

When the above requirements are met, you can compile with:

``
g++ hashmap.cpp -o HMPaint.exe -lz -fopenmp -lpthread -larmadillo
``

To run **HMPaint**, enter the following command:

``
./HMPaint.exe [-parameter1 value1 -parameter2 value2 ......]
``

# Parameters

## Required Parameters

**HMPaint** has below 5 required parameters, all of which are files.

* **-reffile [file]** Reference phase (or gzipped phase) file that contains the genotype data for each reference sample.

* **-targetfile [file]** Target phase (or gzipped phase) file that contains the genotype data for each target sample.

* **-mapfile [file]** Genetic map file that contains two columns with the first line specifying the column names. The first column is the SNP position (in base) and the second column is the genetic distance of each SNP (in Morgan). The number of SNPs must be the same as that in donorfile and targetfile.

* **-popfile [file]** Population file of reference individuals that contains two columns. The first column is the names of reference samples (must be in the same order as donorfile. The second column is the population indices of the reference samples. The population indices must be non-positive integers ranging from 0 to k-1, assuming there are k different populations in the reference panel.

* **-targetname [file]** Target name file that contains the names of target samples. This parameter is necessary because the phase file doesn't contain the sample names.

## Optional Parameters

* **-out [string]** Prefix of the output file names (**default=HMPaint**).

* **-haploid [1/0]** The individuals are haploid (**1**) or diploid (**0**) (**default=0**).

* **-run [paint/chunklength/both]** Run painting and/or LDAS and AAS (**paint**), chunk length of reference panel (**chunklength**) or both analysis (**both**) (**default=paint**).

* **-painting [1/0]** Output the painting results (probabilities) for each individual at each SNP (**1**) or not (**0**) (**default=1**).

* **-aveSNPpainting [1/0]** Output the average painting probabilities for each SNP (**1**) or not (**0**) (**default=1**).

* **-aveindpainting [1/0]** Output the average painting probabilities for each individual (**1**) or not (**0**) (**default=1**).

* **-LDA [1/0]** (**default=1**) Output the LDA results (**1**) or not (**0**) (**default=0**). It might be slow: the computational time is proportional to the number of reference populations and the density of SNPs in the chromosome.

* **-LDAS [1/0]** (**default=1**) Output the LDAS results (**1**) or not (**0**) (**default=0**). It might be slow: the computational time is proportional to the number of reference populations and the density of SNPs in the genome.

* **-AAS [1/0]** (**default=1**) Output the AAS results (**1**) or not (**0**) (**default=1**).

* **-ncores [integer&ge;0]** The number of CPU cores used for the analysis (**default=0**). The default **ncores** parameter uses all the available CPU cores of your device.

* **-L_initial [integer>0]** The initial length of matches (the number of SNPs) that **HMPaint** searches for (**default=320**). **L_initial** must be bigger than **L_minmatch** and should be a power of 2 of **L_minmatch** for computational efficiency.

* **-matchfrac [number&isin;(0,1)]** The proportion of matches of at least **L_minmatch** SNPs that **HMPaint** searches for (**default=0.002**). Positions with more than **matchfrac** proportion of matches of at least **L_minmatch** SNPs will retain at least the longest **matchfrac** proportion of matches. A larger **matchfrac** increases both the accuracy and the computational time.

* **-L_minmatch [integer>0]** The minimal length of matches that **HMPaint** searches for (**default=40**). Positions with fewer than **matchfrac** proportion of matches of at least **L_minmatch** SNPs will retain all the matches of at least **L_minmatch**. A larger **L_minmatch** increases both the accuracy and the computational time.

* **-method [Viterbi/EM]** The algorithm used for estimating the recombination scaling constant (**default=Viterbi**).

* **-diff_rho [1/0]** Use different recombination scaling constant (**1**) or the same value (**0**) for each target sample (**default=0**).

* **-fixrho [number&ge;0]** The value of the fixed recombination scaling constant (**default=0**). **HMPaint** will estimate rho as the average recombination scaling constant of **indfrac** target samples under the default **fixrho** and **diff_rho**.

* **-indfrac [number&isin;(0,1)]** The proportion of individuals used to estimate the recombination scaling constant (**default=0.1**).

* **-minsnpEM [integer>0]** The minimum number of SNPs used for EM algorithm if **-method EM** is specified (**default=2000**).

* **-EMsnpfrac [number&isin;(0,1)]** The proportion of SNPs used for EM algorithm if **-method EM** is specified (**default=0.1**). Note that if nsnp***EMsnpfrac** < **minsnpEM**, **minsnpEM** SNPs will be used for EM algorithm.

* **-ite_time [integer>0]** (**default=10**)

* **-window [number>0]** (**default=0.04**)

* **-LDAfactor [integer&ge;1]** (**default=1**)
