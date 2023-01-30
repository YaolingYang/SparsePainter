# To make sure files contain the same SNPs, we need to run
vcftools --gzvcf p2.vcf.gz --positions-overlap p1.vcf.gz --out p2new --recode
vcftools --vcf p2new.recode.vcf --positions-overlap p3.vcf.gz --out p2new --recode

vcftools --gzvcf p1.vcf.gz --positions-overlap p2new.recode.vcf --out p1new --recode
vcftools --gzvcf p3.vcf.gz --positions-overlap p2new.recode.vcf --out p3new --recode

for i in `seq 1 3`; do mv p"$i"new.recode.vcf p"$i"new.vcf; done

# Longmatchquery -- find matches longer than L.
# take p1new.vcf as panel haplotypes, and p3new.vcf as query haplotypes.

for i in {1..2}; do ./exelmq_swp_dpbwt.exe -i p"$i"new.vcf -q p3new.vcf -m -L 1 -o p"$i"match.txt; done

# we can also use pbwt in longmatchquery, which is faster:
for i in {1..2}; do ./exelmq_swp_pbwt.exe -i p"$i"new.vcf -q p3new.vcf -m -L 1 -o p"$i"match.txt; done


# we want to compare the speed of d-pbwt software with pbwt software
# the pbwt software is only able to find the set-maximal matches

pbwt -readVcfGT p1new.vcf -write p1.pbwt
pbwt -readVcfGT p2new.vcf -write p2.pbwt
pbwt -readVcfGT p3new.vcf -write p3.pbwt

pbwt -read p1.pbwt -matchIndexed p3.pbwt > p1matchp3.txt
pbwt -read p2.pbwt -matchIndexed p3.pbwt > p2matchp3.txt
