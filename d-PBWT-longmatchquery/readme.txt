This program takes VCF files as input. It outputs matches between query haplotypes and haplotypes in the panel length L or longer. It is recommended for multiallelic sites to be removed from the input VCF. If multiallelic sites are not removed, the code will take 0 as 0 and any other allele as 1. The match file output is of the form "X = qY at [a,b)". This means haplotype X in the input vcf matches query haplotype Y in the shuffled order from site a to site b, b excluded. Site 0 is the top most site in the vcf and haplotype 0 is the left most haplotype in the vcf. The algorithms used are presented in https://www.biorxiv.org/content/10.1101/2020.01.14.906487v1.

Compile with std=c++17 or higher. The following command may be used:
g++ -O3 -std=c++17 -o exelmq_xswp_xpbwt.exe lmq_xswp_xpbwt.cpp

Where xswp is 3swp or 1swp and xpbwt is dpbwt or pbwt.

Sample Usage:

After compilation, you can run the program on the sample vcfs using the following command:
./exelmq_xwp_xpbwt.exe -i example.vcf -q query_example.vcf -m -L 3

To generate a random order and only use only input file with the last two haplotypes as query haplotypes:
./exelmq_xwp_xpbwt.exe -i example.vcf -n 2 -L 3 -g order.txt

To use the previously generated order with the last 5 sequences in the shuffled order as query haplotypes:
./exelmq_xwp_xpbwt.exe -i example.vcf -n 5 -L 3 -r order.txt 

options:
-i, input VCF file for the panel, default is panel.vcf
-m, use separate input files for panel and query haplotypes (default off), if not passed, -n must be specified
-n, number of query haplotypes (only if query haplotypes are in the input vcf, the last n haplotypes in the vcf will be used), default is 200
-o, output file for matches, default is matches.txt
-t, output file for time, default is longMatchTime.txt
-L, length of match to output. Matches length L or longer will be outputted, default is 1000
-q, input VCF for query haplotypes to search against the panel, default is query.vcf
-r, read and input order from file provided.
-g, generate and use a randomized order, a file name must be provided with this option to write the order used  
-d, use for debugging, this outputs the haplotype, prefix, divergence, u, and v panels. This is meant to only be used with small panels
-h, help
-z, run with no parameters
no parameters, help
