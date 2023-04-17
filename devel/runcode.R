library("Rcpp")
Sys.setenv(PKG_CXXFLAGS = "-I.")
Sys.setenv(PKG_LIBS = paste("-lz", paste0(normalizePath(getwd()), "/gzstream.o"), sep = " "))
sourceCpp("hashmap.cpp", env = globalenv(), rebuild = TRUE, verbose = TRUE)

setwd("data")

a=paintingalldense(nind=50,targetfrac=1,fixrho=FALSE,method="Viterbi",
                   targetfile="target.vcf.gz",donorfile="donor.vcf.gz",
                   minmatchfrac = 0.002,L_initial=320,L_minmatch=20,window=0.05)
