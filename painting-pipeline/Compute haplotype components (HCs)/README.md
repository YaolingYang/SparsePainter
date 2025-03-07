# Pipeline to compute haplotype components (HCs)
Here we describe the pipeline for computing haplotype components (HCs) through [PBWTpaint](https://github.com/richarddurbin/pbwt). We should have a genotype file from each chromosome, i.e. 22 genotype files in total. Let i denote the chromosome index, the input files are named as ``chr${i}_UKBall.vcf.gz``.

## Step 1: Run PBWTpaint for each chromosome  
The first step is to run the below command for each chromosome (i from 1 to 22, you may split it into 22 array jobs):  

```pbwt -readVcfGT chr${i}_UKBall.vcf.gz -paintSparse chr${i}_UKBall 100 2 500```

The explanation of each parameter can be found by typing `pbwt`. In specific, the last parameter controls the sparsity of matches (the larger, the sparser), which is important for datasets with large numbers of individuals, such as the UK Biobank.

(for other input formats please just follow the pbwt instructions on reading data)  

The last parameter controls the sparsity.  

This command generates ``chr${i}_UKBall.chunklengths.s.out.gz`` for i from 1 to 22.  

## Step 2: Generate the overall chunk length matrix

The next step is to do a "weighted" sum of each entry of the N*N sparse matrix (N is the number of individuals), because PBWTpaint reports the chunk length based on the number of SNPs while it should be based on the genetic distance. Here we provide an example **C++** code ``combine_chunklength.cpp``. Please **replace row 197-198 with the number of SNPs for each chromosome in the dataset**, (row 199-201 is the genetic distance in centiMorgan for each chromosome, which do not need to change if you use a standard recombination map), and **in line 203 replace 487409 with the actual N**. If you use a different file name, please also **replace row 214 and 218 with the actual file name**.    

Note that to run this step, we require chunklength files ``chr${i}_UKBall.chunklengths.s.out.gz`` to be generated from ``-paintSparse`` command form ``pbwt``, which have 3 columns: the first two columns are individual index (integer, starting from 1), and the 3rd column is the chunk length.  

Please ensure ``gzstream.C`` and ``gzstream.h`` are in the same directory as ``combine_chunklength.cpp``, and then compile with:

```g++ combine_chunklength.cpp -o combine -lz -lpthread -llapack -lblas -std=c++0x -g -O3```

After which you run with

```./combine```

This may take few hours to run, and afterwards we obtain the output ``full_chunklength_UKBall.txt.gz``. The first two columns are the individual indices (IND1 and IND2), and the third column is the expected chunk length that IND1 is copied from IND2 in centiMorgan. This output file represents the sparse chunk length matrix from all-vs-all painting.


## Step 3: Compute HCs via SVD

The last step is to do SVD to obtain HCs. Here we provide an example code in R, please **replace nsnp=487409 with the actual N**. 

Note that if ``full_chunklength_UKBall.txt.gz`` has more than 2^31 rows (the limit of R), then we need to remove some weakly associated individual pairs (i.e. removing the rows with the smallest numbers of the last column of ``full_chunklength_UKBall.txt.gz``) and re-weight (each row of the sparse matrix should sum up to be the total genetic distance in centiMorgan, which is 3545.04 in the standard genetic map that we use), with e.g. data chunking approach.

```
library(data.table)  
library(Matrix)  
library(sparsesvd)  
nsnp=487409  
number_of_HCs=100  
cat("Begin reading data\n")  
cl <- fread("full_chunklength_UKBall.txt.gz")  
cat("Begin making sparse matrix\n")  
A <- sparseMatrix(cl$V1, cl$V2, x = log10(cl$V3+1), dims=c(nsnp,nsnp))  
cat("Begin svd \n")  
res<-sparsesvd(A,rank=number_of_HCs)  
cat("Begin calculate HCs\n")  
HCs <- res$u %*% diag(sqrt(res$d))  
cat("Begin writing HCs\n")  
colnames(HCs)=paste0(“HC”,1:number_of_HCs)  
fwrite(HCs, "HCs_UKBall.csv", row.names = FALSE, sep = ',')  
```

Now we get the top 100 HCs in file ``HCs_UKBall.csv``. The individual orders are the same as the original vcf.gz file.

Please also be aware that there should exist more efficient ways to do Step 2 and 3.
