# HMPaint
This is the location for painting with d-PBWT.

The main code is in devel/hashmap.cpp.

To run the code, you should load the [Armadillo](https://arma.sourceforge.net/download.html) and [boost](https://www.boost.org/) libraries, and also have ["gzstream.h" and "gzstream.C"](https://www.cs.unc.edu/Research/compgeom/gzstream/) in your directory. 

When the above requirements are met, you can compile with:

``
g++ hashmap.cpp -o test.exe -lz -fopenmp -lpthread -L/mnt/storage/software/libraries/gnu/12.4.0/lib64 -larmadillo
``

Then, run the codes with (for example)

``
chr=6
pop='british'
./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_${pop}.phase.gz -targetname ../${pop}.txt -avepaintingfile avepainting_chr${chr}_${pop}.txt.gz -LDAfile LDA_chr${chr}_${pop}.txt.gz -LDASfile LDAS_chr${chr}_${pop}.txt -AASfile AAS_chr${chr}_${pop}.txt -mapfile chr${chr}map.txt -popfile ../popnames.txt
``
