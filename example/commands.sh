#!/bin/bash
# To make sure files contain the same SNPs, we need to run
exe=""
bindir="../d-PBWT-longmatchquery/"
if [ "$OS" = "Windows_NT" ]; then
    detected_OS="Windows"
    exe=".exe"
    bindir="..\\d-PBWT-longmatchquery\\"
else
    detected_OS=`uname 2>/dev/null || echo Unknown`
fi
echo "Using OS $detected_OS"

## 
vcftools --gzvcf p2.vcf.gz --positions-overlap p1.vcf.gz --out p2new --recode
mv p2new.recode.vcf p2new.recode.tmp.vcf ## Race condition: can't write to input file in vcftools
vcftools --vcf p2new.recode.tmp.vcf --positions-overlap p3.vcf.gz --out p2new --recode

vcftools --gzvcf p1.vcf.gz --positions-overlap p2new.recode.vcf --out p1new --recode
vcftools --gzvcf p3.vcf.gz --positions-overlap p2new.recode.vcf --out p3new --recode

for i in `seq 1 3`; do mv p${i}new.recode.vcf p${i}new.vcf; done

for i in {1..3}; do bgzip -k -f p${i}new.vcf; tabix p${i}new.vcf.gz; done

for i in {1..3}; do bcftools norm -d all p${i}new.vcf.gz | bgzip -f > p${i}new.bak.vcf.gz; mv p${i}new.bak.vcf.gz p${i}new.vcf.gz; tabix -f p${i}new.vcf.gz; done

for i in {1..3}; do zcat p${i}new.vcf.gz > p${i}new.vcf ; done

# Longmatchquery -- find matches longer than L.
# take p1new.vcf as panel haplotypes, and p3new.vcf as query haplotypes.

for i in {1..2}; do ${bindir}exelmq_1swp_dpbwt$exe -i p${i}new.vcf -q p3new.vcf -m -L 1 -o p${i}match.txt; done

# we can also use pbwt in longmatchquery, which is faster but results may be incorrect!
# for i in {1..2}; do ${bindir}exelmq_1swp_pbwt$exe -i p${i}new.vcf -q p3new.vcf -m -L 1 -o p${i}match.txt; done

# Then we merge the match files into a single file, which is the input to R
for i in {2..2}; 
do awk -v var="$(((i-1)*40))" '{$1=$1+var; print $0}' p${i}match.txt > p${i}matchnew.txt; 
   mv p${i}matchnew.txt p${i}match.txt; 
   cat p$((i-1))match.txt p${i}match.txt > p${i}matchnew.txt;
   mv p${i}matchnew.txt p${i}match.txt;
done

mv p2match.txt matchdata.txt

# we need the map showing the recombination distance
plink --vcf p1new.vcf --recode --out p1new

# split it into a match file for a target individual
#for i in {0..999}; do grep -w "q$i" matchdata.txt > match$i.txt; done

# we can use parallel running:
seq 0 999 | parallel --jobs 10 "grep -w q{} matchdata.txt > match{}.txt"


# we want to compare the speed of d-pbwt software with pbwt software
# the pbwt software is only able to find the set-maximal matches

#pbwt -readVcfGT p1new.vcf.gz -write p1.pbwt
#pbwt -readVcfGT p2new.vcf.gz -write p2.pbwt
#pbwt -readVcfGT p3new.vcf.gz -write p3.pbwt

#pbwt -read p1.pbwt -matchIndexed p3.pbwt > p1matchp3.txt
#pbwt -read p2.pbwt -matchIndexed p3.pbwt > p2matchp3.txt
