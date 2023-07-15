# HMPaint
This is the location for painting with d-PBWT.

The main code is in devel/hashmap.cpp.

To run the code, you should load the [Armadillo](https://arma.sourceforge.net/download.html) library, and also have ["gzstream.h" and "gzstream.C"](https://www.cs.unc.edu/Research/compgeom/gzstream/) in your directory. 

When the above requirements are met, you can compile with:

``
g++ hashmap.cpp -o test.exe -lz -fopenmp -lpthread -L/mnt/storage/software/libraries/gnu/12.4.0/lib64 -larmadillo
``

# Parameters

