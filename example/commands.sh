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
bcftools index --force p1.vcf.gz
bcftools index --force p2.vcf.gz
bcftools index --force p3.vcf.gz

# Use bcftools isec to find overlapping SNPs
bcftools isec -n=3 p1.vcf.gz p2.vcf.gz p3.vcf.gz -o common.vcf

# Use bcftools view to extract overlapping SNPs from the original .vcf.gz files
bcftools view -R common.vcf p1.vcf.gz > p1_common.vcf
bcftools view -R common.vcf p2.vcf.gz > p2_common.vcf
bcftools view -R common.vcf p3.vcf.gz > p3_common.vcf

for i in {1..3}; do bgzip -k -f p${i}_common.vcf; tabix p${i}_common.vcf.gz; done

for i in {1..3}; do bcftools norm -d all p${i}_common.vcf.gz | bgzip -f > p${i}_common.bak.vcf.gz; mv p${i}_common.bak.vcf.gz p${i}_common.vcf.gz; tabix -f p${i}_common.vcf.gz; done

for i in {1..3}; do zcat p${i}_common.vcf.gz > p${i}_common.vcf ; done

bcftools merge p1_common.vcf.gz p2_common.vcf.gz  --force-samples -o p_donor.vcf

mv p3_common.vcf p_target.vcf

for i in `seq 1 40`; do echo "p1:$i p1 1" >> popnames.ids; done
for i in `seq 1 40`; do echo "p2:$i p2 1" >> popnames.ids; done
for i in `seq 1 1000`; do echo "$i p_target 1" >> popnames.ids; done

# The above generates the required input to the software
# The below begins running the software, with the input data of p_donor.vcf, p_target.vcf and popnames.ids

# Longmatchquery -- find matches longer than L.
${bindir}./exelmq_1swp_dpbwt$exe -i p_donor.vcf -q p_target.vcf -m -L 1 -o target_match.txt
${bindir}./exelmq_1swp_dpbwt$exe -i p_donor.vcf -q p_donor.vcf -m -L 1 -o donor_match.txt

# we need the map showing the recombination distance
plink --vcf p_donor.vcf --recode --out p

# split it into a match file for a target/donor individual
seq 0 999 | parallel --jobs 10 "grep -w q{} target_match.txt > target_match{}.txt"
seq 0 79 | parallel --jobs 10 "grep -w q{} donor_match.txt > donor_match{}.txt"
