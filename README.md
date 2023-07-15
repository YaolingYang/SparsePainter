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

**HMPaint** has below 5 required parameters.

* **-donorfile <file>** reference phase (or gzipped phase) file that contains the genotype data for each reference (donor) sample.

* **-targetfile <file>** target phase (or gzipped phase) file that contains the genotype data for each target sample.

* **-mapfile <file>** recombination file which contains two columns with the first line specifying the column names. The first column is the SNP position (in base) and the second column is the genetic distance of each SNP (in Morgan). The number of SNPs must be the same as that in donorfile and targetfile.

*  **-popfile <file>** 

