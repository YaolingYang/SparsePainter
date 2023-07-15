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

* **-donorfile [file]** Reference phase (or gzipped phase) file that contains the genotype data for each reference (donor) sample.

* **-targetfile [file]** Target phase (or gzipped phase) file that contains the genotype data for each target sample.

* **-mapfile [file]** Genetic map file that contains two columns with the first line specifying the column names. The first column is the SNP position (in base) and the second column is the genetic distance of each SNP (in Morgan). The number of SNPs must be the same as that in donorfile and targetfile.

* **-popfile [file]** Population file of reference individuals that contains two columns. The first column is the names of reference samples (must be in the same order as donorfile. The second column is the population indices of the reference samples. The population indices must be non-positive integers ranging from 0 to k-1, assuming there are k different populations in the reference panel.

* **-targetname [file]** Target name file that contains the names of target samples. This parameter is necessary because the phase file doesn't contain the sample names.

## Optional Parameters

* **-out [string]** Prefix of the output file names (**default=HMPaint**).

* **-run [string]** (**default=paint**)

* **-outputpainting [1/0]** (**default=1**)

* **-outputaveSNPpainting [1/0]** (**default=1**)

* **-outputaveindpainting [1/0]** (**default=1**)

* **-outputLDA [1/0]** (**default=1**)

* **-outputLDAS [1/0]** (**default=1**)

* **-outputAAS [1/0]** (**default=1**)

* **-ncores [integer&ge;0]** (**default=0**)

* **-haploid [1/0]** (**default=0**)

* **-method [string]** (**default=Viterbi**)

* **-diff_rho [1/0]** (**default=0**)

* **-fixrho [number&ge0]** (**default=0**)

* **-L_initial [integer>0]** (**default=320**)

* **-minmatchfrac [number&isin;(0,1)]** (**default=0.002**)

* **-L_minmatch [integer>0]** (**default=40**)

* **-indfrac [number&isin;(0,1)]** (**default=0.1**)

* **-minsnpEM [integer>0]** (**default=10000**)

* **-EMsnpfrac [number&isin;(0,1)]** (**default=0.1**)

* **-ite_time [integer>0]** (**default=10**)

* **-window [number>0]** (**default=0.04**)

* **-LDAfactor [integer&ge;1]** (**default=1**)
