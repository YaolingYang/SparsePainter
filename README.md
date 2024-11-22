# SparsePainter

## Note: This experimental fork removes the longest match from consideration in determining the sparse state space of the Li and Stephens HMM.

----------------------

**SparsePainter** is an efficient tool for local ancestry inference (LAI) coded in C++. It extends [**PBWT**](https://github.com/richarddurbin/pbwt) algorithm to find K longest matches at each position, and uses the **Hash Map** structure to implement the forward and backward algorithm in the Hidden Markov Model (HMM) leveraging the sparsity of haplotype matches. SparsePainter can infer **fine-scale local ancestry** (per individual per SNP) and **genome-wide total ancestry**, it also enables efficiently calculating [**Linkage Disequilibrium of Ancestry (LDA), LDA score (LDAS)**](https://github.com/YaolingYang/LDAandLDAscore) and [**Ancestry Anomaly Score (AAS)**](https://github.com/danjlawson/ms_paper) for understanding the population structure, evolution, selection, etc..  

-   Authors:  
    Yaoling Yang (<yaoling.yang@bristol.ac.uk>)  
    Daniel Lawson (<dan.lawson@bristol.ac.uk>)

-   Maintainer:  
    Yaoling Yang (<yaoling.yang@bristol.ac.uk>)  

-   **SparsePainter website:**  https://sparsepainter.github.io/

-   Version: 1.2.1 (**[Changelog](#changelog)**)

-   **SparsePainter and PBWTpaint Reference:** [Yang, Y., Durbin, R., Iversen, A.K.N & Lawson, D.J. Sparse haplotype-based fine-scale local ancestry inference at scale reveals recent selection on immune responses. medRxiv (2024). doi:10.1101/2024.03.13.24304206](https://www.medrxiv.org/content/10.1101/2024.03.13.24304206v2)

-   **Pipeline for** [biobank-scale painting](https://github.com/YaolingYang/SparsePainter/tree/main/painting-pipeline/standard%20painting) and [computing haplotype components (HCs)](https://github.com/YaolingYang/SparsePainter/tree/main/painting-pipeline/Compute%20haplotype%20components%20(HCs)) **are available.**

-   **LDA, LDA score and AAS Reference:** [Barrie, W., Yang, Y., Irving-Pease, E.K. et al. Elevated genetic risk for multiple sclerosis emerged in steppe pastoralist populations. Nature 625, 321–328 (2024)](https://www.nature.com/articles/s41586-023-06618-z)

-   **PBWTpaint GitHub Repository:** [https://github.com/richarddurbin/pbwt](https://github.com/richarddurbin/pbwt)

-   Overview of SparsePainter and PBWTpaint
![overview](overview.png)

# Installation

To install SparsePainter, please follow the below steps.  
``git clone git@github.com:YaolingYang/SparsePainter.git``  
``cd SparsePainter``  
``make``  

To update the newer version of SparsePainter, you can remove lines 10-12 of Makefile, since armadillo has already been installed during your initial installation.

# Dependencies

SparsePainter requires g++ >=6 and depends on   
[Armadillo-v12.6.5](https://arma.sourceforge.net/download.html) to compute AAS;   
[gzstream-v1.5](https://www.cs.unc.edu/Research/compgeom/gzstream/) to read and write gzipped files.


# Usage

Either variant call format (VCF) or phase format is supported by **SparsePainter**. Both files should be phased and without missing data. Inputting phase format is slightly faster than inputting the VCF format. To prepare the phase format for **SparsePainter**, you should get [PBWT](https://github.com/richarddurbin/pbwt) installed, which converts Variant Call Format (VCF) to phase format by the following command:

``
pbwt -readVcfGT XXX.vcf -writePhase XXX.phase
``

To run **SparsePainter**, enter the following command:

``
./SparsePainter [-command1 -command2 ...... -command3 parameter3 -command4 parameter4 ......]
``


# Commands

## Required Commands

**SparsePainter** has below 6 required commands together with additional commands that specify the desired output.

* **-reffile [file]** Reference vcf (including gzipped vcf), or phase (including gzipped phase) file that contains the (phased non-missing) genotype data for the reference samples.

* **-targetfile [file]** Reference vcf (including gzipped vcf), or phase (including gzipped phase) file that contains the (phased non-missing) genotype data for the target samples. To paint reference samples against themselves, please set ``targetfile`` to be the same as ``reffile``. The file type of ``targetfile`` and ``reffile`` should be the same.

* **-mapfile [file]** Genetic map file that contains two columns **with headers**. The first column is the SNP position (in base) and the second column is the genetic distance of each SNP (in **centiMorgan**). The SNPs must be the same and of the same order as those in ``reffile`` and ``targetfile``.

* **-popfile [file]** Population file of reference individuals that contains two columns **without headers**. The first column is the names of all the reference samples (must be in the same order as ``reffile``). The second column is the population labels of the reference samples, which can be either strings or numbers.

* **-namefile [file]** Name file that contains the names of samples to be painted, following the same order as they appear in ``targetfile``.

* **-out [string]** Prefix of the output file names (**default=SparsePainter**).

**At least one of the below commands should also be given in order to run SparsePainter**

* **-prob** Output the local ancestry probabilities for each target sample at each SNP. The output is a gzipped text file (.txt.gz) with format specified in `-probstore`.

* **-chunklength** Output the expected length (in centiMorgan) of copied chunks of each local ancestry for each target sample. The output is a gzipped text file (.txt.gz).

* **-chunkcount** Output the expected number of copied chunks of each local ancestry for each target sample. The output is a gzipped text file (.txt.gz).

* **-aveSNP** Output the average local ancestry probabilities for each SNP. The output is a text file (.txt).

* **-aveind** Output the average local ancestry probabilities for each target individual. The output is a text file (.txt).

* **-LDA** Output the Linakage Disequilibrium of Ancestry (LDA) of each pair of SNPs. The output is a gzipped text file (.txt.gz). It might be slow: the computational time is proportional to the number of local ancestries and the density of SNPs in the chromosome.

* **-LDAS** Output the Linakage Disequilibrium of Ancestry Score (LDAS) of each SNP. The output is a text file (.txt), including the LDAS and its lower and upper bound, which can be used for quality control. It might be slow: the computational time is proportional to the number of local ancestries and the density of SNPs in the genome.

* **-AAS** Output the test statistic of Ancestry Anomaly Score (AAS) of each SNP. The output is a text file (.txt). The AAS test statistic follows chi-squared distribution with K degrees of freedom under the null, where K is the number of reference populations.

## Optional Commands

### Commands without values

* **-haploid** The individuals are haploid.

* **-diff_lambda** Use different recombination scaling constants for each target sample. If this command is not given, the fixed lambda will be output in a text file (.txt) for future reference.

* **-loo** Paint with leave-one-out strategy: one individual is left out of each population (self from own population). If `-loo` is not specified under reference-vs-reference painting (`reffile=targetfile`), each individual will be automatically left out of painting. For accuracy, please do not use this command if any of the reference populations has very few (e.g. <=5) samples.

* **-rmrelative** Leave out the reference sample that is the most related to the target sample under leave-one-out mode (`-loo`), if they share at least ``relafrac`` proportion of SNPs of a continuous segment. Please do not use this command for reference-vs-reference painting.

* **-outmatch** Output the number of matches at each SNP for each target haplotype. The output file format is a gzipped text file (.txt.gz).

### Commands with values

* **-ncores [integer&ge;0]** The number of CPU cores used for the analysis (**default=0**). The default ``ncores`` uses all the available CPU cores of your device.

* **-fixlambda [number&ge;0]** The value of the fixed recombination scaling constant (**default=0**). **SparsePainter** will estimate lambda as the average recombination scaling constant of ``indfrac`` target samples under the default ``fixlambda`` and ``diff_lambda``.

* **-nmatch [integer>=1]** The number of haplotype matches of at least ``Lmin`` SNPs that **SparsePainter** searches for (**default=10**). Positions with more than ``nmatch`` matches of at least ``Lmin`` SNPs will retain at least the longest ``nmatch`` matches. A larger ``nmatch`` slightly improves accuracy but significantly increases the computational time.

* **-L0 [integer>0]** The initial length of matches (the number of SNPs) that **SparsePainter** searches for (**default=320**). ``L0`` must be bigger than ``Lmin`` and preferrably be a power of 2 of ``Lmin`` for computational efficiency.

* **-Lmin [integer>0]** The minimal length of matches that **SparsePainter** searches for (**default=20**). Positions with fewer than ``nmatch`` matches of at least ``Lmin`` SNPs will retain all the matches of at least ``Lmin``. A larger ``Lmin`` increases both the accuracy and the computational time.

* **-method [Viterbi/EM]** The algorithm used for estimating the recombination scaling constant (**default=Viterbi**).

* **-probstore [raw/constant/linear/cluster]** Output the local ancestry probabilities in raw, constant, linear or cluster form (**default=constant**). For each haplotype, in ``raw`` form, we output the probabilities of each SNP with the SNP name being their physical positions in base; in ``constant`` form, we output the range of SNP index, and the painting probabilities that those SNPs share; in ``linear`` form, we output the range of SNP index, and the painting probabilities of the start SNP and the end SNP, while the intermediate SNPs are estimated by the simple linear regression with root mean squared error smaller than ``rmsethre``; in ``cluster`` form, we perform K-means clustering on the painting of each haplotype with ``ncluster`` clusters and maximum ``max_ite`` iterations, and output the average probabilities of each cluster and the cluster of each SNP. Storing in ``constant`` considerably reduces the file size while has the same accuracy compared with storing in ``raw``; storing in ``linear`` has an even smaller file size but becomes slightly slower and loses some accuracy; storing in ``cluster`` has the smallest file size but with slowest speed.

* **-dp [integer>0]** The decimal places of the output of local ancestry probabilities (**default=2**). This also controls the size of the output file for local ancestry probabilities.

* **-rmsethre [number&isin(0,1)]** The upper bound that the root mean squared error of the estimated local ancestry probabilities (**default=0.01**) when storing them in linear form by argument, i.e. ``-probstore linear``.

* **-relafrac [number&isin;(0,1)]** The proportion of the total number of SNPs shared between a reference and target haplotype sample (**default=0.2**). The reference sample will be removed under the leave-one-out (``-loo``) and remove relative (``-rmrelative``) modes.

* **-ncluster [integer>0]** The number of clusters (**default=100**) for K-means clustering under ``-probstore cluster`` mode.

* **-kmeans_ite [integer>0]** The number of maximum iterations (**default=30**) for K-means clustering under ``-probstore cluster`` mode.

* **-SNPfile [file]** File contains the specific physical position (in base) of the SNPs whose local ancestry probabilities are output in the ``raw`` form. If this file is not specified (default), then all the SNPs' local ancestry probabilities will be output in the form specified by ``probstore``. 

* **-indfrac [number&isin;(0,1]]** The proportion of individuals used to estimate the recombination scaling constant (**default=0.1**).

* **-minsnpEM [integer>0]** The minimum number of SNPs used for EM algorithm if ``-method EM`` is specified (**default=2000**).

* **-EMsnpfrac [number&isin;(0,1]]** The proportion of SNPs used for EM algorithm if ``-method EM`` is specified (**default=0.1**). Note that if ``nsnp*EMsnpfrac < minsnpEM``, ``minsnpEM`` SNPs will be used for EM algorithm.

* **-EM_ite [integer>0]** The iteration times for EM algorithm if ``-method EM`` is specified (**default=10**).

* **-window [number>0]** The window for calculating LDA score (LDAS) in centiMorgan (**default=4**).

* **-matchfile [file]** The file name of the set-maximal match file which is the output of [pbwt -maxWithin](https://github.com/richarddurbin/pbwt). This can only be used for painting reference samples against themselves. When ``matchfile`` is given, there is no need to provide ``reffile`` and ``targetfile``, because all the match information required for painting is contained in ``matchfile``. Using set-maximal matches is not recommended because set-maximal matches are extremely sparse and will significantly reduce the accuracy, despite saving compute time.


# Generate mapfile with fixed recombination rate
When analysing human genomes, we suggest using the real human recombination map (available from various online resources) which has varying recombination rates throughout the genome. Please follow the file formats of `mapfile` for generating it.  

In certain scenarios, for example when analysing non-human species, you may use a fixed recombination rate throughout the genome. Below we provide the commands for generating the ``mapfile`` from ``input.vcf`` with a fixed recombination rate of 1e-6 cM/bp.  

```
bcftools query -f '%POS\n' input.vcf > sites.txt  
echo -e "pd\tgd" > map.txt  
awk '{print $1"\t"1e-6*$1}' sites.txt >> map.txt  
```

# Examples
Here we provide examples to run SparsePainter. Examples are explained in more detail on our [SparsePainter website](https://sparsepainter.github.io/).  

The example dataset is contained in the /example folder. This example includes 8000 reference individuals from 4 populations with 2091 SNPs (Both vcf version ``donor.vcf.gz`` and phase version ``donor.phase.gz`` are available), and the aim is to paint 500 target individuals (Both vcf version ``target.vcf.gz`` and phase version ``target.phase.gz`` are available). Remember we have compiled SparsePainter in ``SparsePainter``, then we can paint with the following command:

*  If your input file is in vcf or vcf.gz format:

``
./SparsePainter -reffile donor.vcf.gz -targetfile target.vcf.gz -popfile popnames.txt -mapfile map.txt -namefile targetname.txt -out target_vs_ref -prob -chunklength -chunkcount -aveSNP -aveind
``

*  If your input file is in phase of phase.gz format:

``
./SparsePainter -reffile donor.phase.gz -targetfile target.phase.gz -popfile popnames.txt -mapfile map.txt -namefile targetname.txt -out target_vs_ref -prob -chunklength -chunkcount -aveSNP -aveind
``

The output file for this example includes ``target_vs_ref_prob.txt.gz``, ``target_vs_ref_chunklength.txt.gz``, ``target_vs_ref_chunkcount.txt.gz``, ``target_vs_ref_aveSNPprob.txt``, ``target_vs_ref_aveindprob.txt`` and ``target_vs_ref_fixedlambda.txt``.

To paint the reference individuals against themselves with leave-one-out strategy, run with:

*  If your input file is in vcf or vcf.gz format:

``
./SparsePainter -reffile donor.vcf.gz -targetfile donor.vcf.gz -popfile popnames.txt -mapfile map.txt -namefile refname.txt -out ref_vs_ref -prob -chunklength -chunkcount -aveSNP -aveind -loo
``

*  If your input file is in phase or phase.gz format:

``
./SparsePainter -reffile donor.phase.gz -targetfile donor.phase.gz -popfile popnames.txt -mapfile map.txt -namefile refname.txt -out ref_vs_ref -prob -chunklength -chunkcount -aveSNP -aveind -loo
``

The output file for this example includes ``ref_vs_ref_prob.txt.gz``, ``ref_vs_ref_chunklength.txt.gz``, ``ref_vs_ref_chunkcount.txt.gz``, ``ref_vs_ref_aveSNPprob.txt``, ``ref_vs_ref_aveindprob.txt`` and ``ref_vs_ref_fixedlambda.txt``.


# Extract the local ancestry probabilities of certain SNPs from the output files
We have provided the C++ and R codes (C++ is much faster than R) to extract the local ancestry probabilities of certain SNPs from different output formats (constant, linear or raw) in the folder ``process_output``, which also includes the example data files.

The .cpp files with suffix ``_hap.cpp`` extract probabilities for each haplotype, and those without ``_hap.cpp`` extract probabilities for each diploid individual, i.e. average over two copies.

For example, to extract the local ancestry probabilities of 5000 SNPs (whose indices are contained in ``chr19_GWAS_SNPs.txt``) stored in constant form with C++, we first compile with

``g++ extract_prob_constant.cpp -o extract_constant -lz -std=c++0x -g -O3``

Then we should specify 4 arguments: `npop` (the number of reference populations), `probfile` (local ancestry probabilities file), `SNPfile` (the indices (starting from 1) of the SNPs whose local ancestry probabilities are to be extracted) and `out` (the prefix of the output file).

``./extract_constant -npop 26 -probfile chr19_1000G_constant_10inds_prob.txt.gz -SNPfile chr19_GWAS_SNPs.txt -out chr19_1000G_constant``

Then we got K files (K is the number of populations), and the probabilities for each reference population (for all the individuals and selected SNPs) are extracted in one gzipped text file. Note that the C++ code averages the probabilities of each individual over its 2 haplotypes (if it's diploid).

It works similarly for extracting paintings in the linear form:

``g++ extract_prob_linear.cpp -o extract_linear -lz -std=c++0x -g -O3``

``./extract_linear -npop 26 -probfile chr19_1000G_linear_10inds_prob.txt.gz -SNPfile chr19_GWAS_SNPs.txt -out chr19_1000G_linear``

We do not suggest storing in the raw form, as it's not memory-efficient.

# Combine the output of multiple subfiles if you split the target files to run SparsePainter
To paint large biobanks, it is suggested to split the target samples into multiple (hundreds of) subfiles and paint them separately, which saves memory and computational time. However, people may prefer merge the results of those analyses into a single file for subsequent analysis. Here we explain how to do this for each output.

* ``-prob`` (for any storage mode): 
Retain the first subfile, and then append the rows (excluding the first two rows) of the other subfiles.  

* ``-chunklength`` and ``-chunkcount``:
Retain the first subfile, and then append the rows (excluding the first row) of the other subfiles. To obtain genome-wide chunk length and chunk count, please sum over all chromosomes. 

* ``-aveindpainting``:
Retain the first subfile, and then append the rows (excluding the first row) of the other subfiles.  To obtain genome-wide average painting for individuals, please compute the weighted average of all chromosomes (weighted by the number of SNPs in each chromosome).

* ``-aveSNPpainting``:
Compute the weighted average of all subfiles (weighted by the number of samples in each subfile).

* ``-LDA`` and ``-LDAS``:
Compute the weighted average of all subfiles (weighted by the number of samples in each subfile).

* ``-AAS``:
AAS cannot be directly merged. To obtain the overall AAS, please run SparsePainter without -AAS, but with -aveSNPpainting. Then compute the weighted average of all subfiles (weighted by the number of samples in each subfile). Then compile ``doAAS.cpp`` (contained in the folder ``process_output``) with below or similar commands (depending on your device, see ``Makefile``):

``g++ -I./armadillo-12.6.5/include doAAS.cpp -o doAAS -lz -fopenmp -lpthread -L./armadillo-12.6.5 -larmadillo -llapack -lblas -std=c++0x -g -O3 -Wl,-rpath=./armadillo-12.6.5``

Finally run the code:

``./doAAS -aveSNPfile [your weighted average aveSNPpainting file] -out [your output file prefix]``

<a id="changelog"></a>
# Changelog
* **2024-10-07 Version 1.2.1**  
Check if the VCF data input includes genotypes that are not 0 or 1, then return the error.
* **2024-06-18 Version 1.2.0**  
Enable string population labels (the 2nd column of ``popfile``).  
In the 2nd column of ``mapfile``, change the unit of genetic distance from Morgan to **centiMorgan**.
* **2024-04-24 Version 1.1.0**  
Enable the output of expected number of copied chunks by command ``-chunkcount``.
* **2024-03-12 Version 1.0.0**  
Release SparsePainter and preprint.
