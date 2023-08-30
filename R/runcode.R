library("Rcpp")
Sys.setenv(PKG_CXXFLAGS = "-I.")
Sys.setenv(PKG_LIBS = paste("-lz", paste0(normalizePath(getwd()), "/gzstream.o"), sep = " "))
sourceCpp("revised_hashmap.cpp", env = globalenv(), rebuild = TRUE, verbose = TRUE)
setwd("bigeg2")

paintingalldense(fixrho=FALSE,method="Viterbi",minmatchfrac = 0.001,
                 L_initial=320,L_minmatch=40,window=0.05,targetname="targetname2.txt",
                 donorfile="donor.phase.gz",targetfile="donor.phase.gz",ncores=12)
#setwd("UKB")
#paintingalldense(fixrho=FALSE,method="Viterbi",minmatchfrac = 0.002,
#                 L_initial=320,L_minmatch=20,window=0.05,targetname="targetname.txt",
#                 mapfile="chr19_map.txt",
#                 donorfile="chr19.phase.gz",targetfile="chr19.phase.gz")

#paintingalldense(fixrho=FALSE,method="Viterbi",minmatchfrac = 0.002,
#                 L_initial=320,L_minmatch=20,window=0.05,targetname="targetname2.txt",
#                 donorfile="donor.phase.gz",targetfile="donor.phase.gz")

#paintingalldense(fixrho=FALSE,method="Viterbi",minmatchfrac = 0.002,
#                 L_initial=320,L_minmatch=20,window=0.05,ncores=5)
#chunklengthall(donorfile="donor.phase.gz")
#paintingalldense(fixrho=FALSE,method="Viterbi",targetfile="target.vcf.gz",
#                           donorfile="donor.vcf.gz",
#                           minmatchfrac = 0.002,L_initial=320,L_minmatch=20,window=0.05)
#chunklengthall(donorfile="donor.phase.gz")

#setwd("UKB")
#paintingalldense(fixrho=FALSE,method="Viterbi",targetfile="chr19_filtered.vcf.gz",
#                 donorfile="chr19_filtered.vcf.gz",mapfile="chr19_map.txt",
#                 popfile="popnames.txt",targetname="targetname.txt",
#                 minmatchfrac = 0.0002,L_initial=560,L_minmatch=35,window=0.05)
