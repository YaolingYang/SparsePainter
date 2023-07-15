# HMPaint
HMPaint is an efficient tool for local ancestry inference (LAI). It uses d-PBWT algorithm to find K longest matches at each position, and uses the Hash Map strategy to implement the forward and backward algorithm in the Hidden Markov Model (HMM) because of the sparsity of haplotype matches. HMPaint incorporates the function for efficiently calculating [Linkage Disequilibrium of Ancestry (LDA), LDA score (LDAS)](https://github.com/YaolingYang/LDAandLDAscore) and [Ancestry Anomaly Score (AAS)](https://github.com/danjlawson/ms_paper) for understanding the population structure and evolution.

# Installation

The main code is in devel/hashmap.cpp.

To run the code, you should load the [Armadillo](https://arma.sourceforge.net/download.html) library, and also have ["gzstream.h" and "gzstream.C"](https://www.cs.unc.edu/Research/compgeom/gzstream/) in your directory. 

When the above requirements are met, you can compile with:

``
g++ hashmap.cpp -o test.exe -lz -fopenmp -lpthread -larmadillo
``

# Parameters

The required
