# SparsePainter
**SparsePainter** is an efficient tool for local ancestry inference (LAI) coded in C++. It improves **PBWT** algorithm to find K longest matches at each position, and uses the **Hash Map** strategy to implement the forward and backward algorithm in the Hidden Markov Model (HMM) because of the sparsity of haplotype matches. SparsePainter incorporates the function for efficiently calculating [Linkage Disequilibrium of Ancestry (LDA), LDA score (LDAS)](https://github.com/YaolingYang/LDAandLDAscore) and [Ancestry Anomaly Score (AAS)](https://github.com/danjlawson/ms_paper) for understanding the population structure, evolution, selection, etc..  

-   Authors:  
    Yaoling Yang (<yaoling.yang@bristol.ac.uk>)  
    Daniel Lawson (<dan.lawson@bristol.ac.uk>)

# Installation

The main code is in **SparsePainter.cpp**.

You should load the [Armadillo](https://arma.sourceforge.net/download.html) library, and also have ["gzstream.h" and "gzstream.C"](https://www.cs.unc.edu/Research/compgeom/gzstream/) in your directory. 

Either converts variant call format (VCF) or phase format is supported by **SparsePainter**. Inputting phase format is slightly faster than inputting the VCF format. To prepare the phase format for **SparsePainter**, you should get [PBWT](https://github.com/richarddurbin/pbwt) installed, which converts Variant Call Format (VCF) to phase format by the following command:

``
pbwt -readVcfGT XXX.vcf -writePhase XXX.phase
``

When the above requirements are met, you can compile with:

``
g++ SparsePainter.cpp -o SparsePainter.exe -lz -fopenmp -lpthread -larmadillo -std=c++0x -g -O3
``

To run **SparsePainter**, enter the following command:

``
./SparsePainter.exe [-parameter1 value1 -parameter2 value2 ......]
``

An example can be found in the **Example** section below.

# Parameters

## Required Parameters

**SparsePainter** has below 6 required parameters.

* **-reffile [file]** Reference vcf (including gzipped vcf), or phase (including gzipped phase) file that contains the genotype data for all the reference samples.

* **-targetfile [file]** Reference vcf (including gzipped vcf), or phase (including gzipped phase) file that contains the genotype data for each target sample. To paint reference samples against themselves, please set ``targetfile`` to be the same as ``reffile``. The file type of ``targetfile`` and ``reffile`` should be the same.

* **-mapfile [file]** Genetic map file that contains two columns with the first line specifying the column names. The first column is the SNP position (in base) and the second column is the genetic distance of each SNP (in Morgan). The number of SNPs must be the same as that in donorfile and targetfile.

* **-popfile [file]** Population file of reference individuals that contains two columns. The first column is the names of reference samples (must be in the same order as ``reffile``). The second column is the population indices of the reference samples. The population indices must be non-positive integers ranging from 0 to k-1, assuming there are k different populations in the reference panel.

* **-targetname [file]** Target name file that contains the names of target samples. This parameter is necessary because the phase file doesn't contain the sample names. When painting reference samples against themselves, this parameter should contain the names of reference samples.

* **-out [string]** Prefix of the output file names (**default=SparsePainter**).

**At least one of the below commands should also be given in order to run SparsePainter**

* **-prob** Output the local ancestry probabilities for each target sample at each SNP. The output file format is a gzipped text file (.txt.gz).

* **-chunklength** Output the chunk length of each local ancestry for each target sample. The output file format is a text file (.txt).

* **-aveSNP** Output the average painting probabilities for each SNP. The output file format is a text file (.txt).

* **-aveind** Output the average painting probabilities for each target individual. The output file format is a text file (.txt).

* **-LDA** Output the LDA of each pair of SNPs. The output file format is a gzipped text file (.txt.gz). It might be slow: the computational time is proportional to the number of local ancestries and the density of SNPs in the chromosome.

* **-LDAS** Output the LDAS of each SNP. The output file format is a text file (.txt). It might be slow: the computational time is proportional to the number of local ancestries and the density of SNPs in the genome.

* **-AAS** Output the AAS of each SNP. The output file format is a text file (.txt).

## Optional Parameters

### Parameters without values

* **-haploid** The individuals are haploid.

* **-diff_lambda** Use different recombination scaling constants for each target sample. If this parameter is not given, the fixed lambda will be output in a text file (.txt) for future reference.\

* **-loo** Paint with leave-one-out strategy: one individual is left out of each population (self from own population). If `-loo` is not specified under reference-vs-reference painting (`reffile=targetfile`), each individual will be automatically left out of painting.

### Parameters with values

* **-ncores [integer&ge;0]** The number of CPU cores used for the analysis (**default=0**). The default **ncores** parameter uses all the available CPU cores of your device.

* **-fixlambda [number&ge;0]** The value of the fixed recombination scaling constant (**default=0**). **SparsePainter** will estimate lambda as the average recombination scaling constant of ``indfrac`` target samples under the default ``fixlambda`` and ``diff_lambda``.

* **-nmatch [integer>=1]** The number of haplotype matches of at least ``L_minmatch`` SNPs that **SparsePainter** searches for (**default=10**). Positions with more than ``nmatch`` matches of at least ``L_minmatch`` SNPs will retain at least the longest ``nmatch`` proportion of matches. A larger ``nmatch`` slightly improves accuracy but significantly increases the computational time.

* **-L0 [integer>0]** The initial length of matches (the number of SNPs) that **SparsePainter** searches for (**default=320**). ``L_initial`` must be bigger than ``L_minmatch`` and should be a power of 2 of ``L_minmatch`` for computational efficiency.

* **-Lmin [integer>0]** The minimal length of matches that **SparsePainter** searches for (**default=20**). Positions with fewer than ``nmatch`` matches of at least ``L_minmatch`` SNPs will retain all the matches of at least ``L_minmatch``. A larger ``L_minmatch`` increases both the accuracy and the computational time.

* **-method [Viterbi/EM]** The algorithm used for estimating the recombination scaling constant (**default=Viterbi**).

* **-indfrac [number&isin;(0,1)]** The proportion of individuals used to estimate the recombination scaling constant (**default=0.1**).

* **-minsnpEM [integer>0]** The minimum number of SNPs used for EM algorithm if ``-method EM`` is specified (**default=2000**).

* **-EMsnpfrac [number&isin;(0,1)]** The proportion of SNPs used for EM algorithm if ``-method EM`` is specified (**default=0.1**). Note that if ``nsnp*EMsnpfrac < minsnpEM``, ``minsnpEM`` SNPs will be used for EM algorithm.

* **-ite_time [integer>0]** The iteration times for EM algorithm if ``-method EM`` is specified (**default=10**).

* **-window [number>0]** The window for calculating LDA score (LDAS) in Morgan (**default=0.04**).

* **-matchfile [file]** The file name of the set-maximal match file which is the output of [pbwt -maxWithin](https://github.com/danjlawson/pbwt/blob/master/pbwtMain.c). This can only be used for painting reference samples against themselves. When ``matchfile`` is given, there is no need to provide ``reffile`` and ``targetfile``, because all the match information required for painting is contained in ``matchfile``. Using set-maximal matches is not recommended because set-maximal matches are extremely sparse that will significantly reduce the accuracy, despite saving compute time.

# Example
The example dataset is contained in the /example folder. This example includes 8000 reference individuals from 4 populations with 2091 SNPs (``donor.phase.gz``), and the aim is to paint 500 target individuals (``target.phase.gz``). Remember we have compiled SparsePainter in ``SparsePainter.exe``, then we can paint with the following command:

(a) If your input file is in vcf or vcf.gz format:

``
./SparsePainter.exe -reffile donor.vcf.gz -targetfile target.vcf.gz -popfile popnames.txt -mapfile map.txt -targetname targetname.txt -out target_vs_ref -prob -chunklength -aveSNP -aveind
``

(b) If your input file is in phase of phase.gz format:

``
./SparsePainter.exe -reffile donor.phase.gz -targetfile target.phase.gz -popfile popnames.txt -mapfile map.txt -targetname targetname.txt -out target_vs_ref -prob -chunklength -aveSNP -aveind
``

The output file for this example includes ``target_vs_ref_prob.txt.gz``, ``target_vs_ref_chunklength.txt.gz``, ``target_vs_ref_aveSNPprob.txt``, ``target_vs_ref_aveindprob.txt`` and ``target_vs_ref_fixedlambda.txt``.

To paint the reference individuals against themselves with leave-one-out strategy, run with:

(a) If your input file is in vcf or vcf.gz format:

``
./SparsePainter.exe -reffile donor.vcf.gz -targetfile donor.vcf.gz -popfile popnames.txt -mapfile map.txt -targetname refname.txt -out ref_vs_ref -prob -chunklength -aveSNP -aveind -loo
``

(b) If your input file is in phase or phase.gz format:

``
./SparsePainter.exe -reffile donor.phase.gz -targetfile donor.phase.gz -popfile popnames.txt -mapfile map.txt -targetname refname.txt -out ref_vs_ref -prob -chunklength -aveSNP -aveind -loo
``

The output file for this example includes ``ref_vs_ref_prob.txt.gz``, ``ref_vs_ref_chunklength.txt.gz``, ``ref_vs_ref_aveSNPprob.txt``, ``ref_vs_ref_aveindprob.txt`` and ``ref_vs_ref_fixedlambda.txt``.
