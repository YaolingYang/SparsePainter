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


exe=""
bindir="../d-PBWT-longmatchquery/"
awk="mawk" # awk is more available but slower
data="data/"
matches="matches/"

if [ "$OS" = "Windows_NT" ]; then
    detected_OS="Windows"
    exe=".exe"
    bindir="..\\d-PBWT-longmatchquery\\"
else
    detected_OS=`uname 2>/dev/null || echo Unknown`
fi
echo "Using OS $detected_OS"
plink="plink1.9"
mkdir -p $matches
mkdir -p $data

## 
for i in {1..3}; do  bcftools index --force p${i}.vcf.gz; done

# Use bcftools isec to find overlapping SNPs
bcftools isec -n=3 p1.vcf.gz p2.vcf.gz p3.vcf.gz -o ${data}common.vcf

# Use bcftools view to extract overlapping SNPs from the original .vcf.gz files
for i in {1..3}; do bcftools view -R ${data}common.vcf p${i}.vcf.gz > ${data}p${i}_common.vcf; done

for i in {1..3}; do bgzip -k -f ${data}p${i}_common.vcf; tabix ${data}p${i}_common.vcf.gz; done

for i in {1..3}; do bcftools norm -d all ${data}p${i}_common.vcf.gz | bgzip -f > ${data}p${i}_common.bak.vcf.gz; mv ${data}p${i}_common.bak.vcf.gz ${data}p${i}_common.vcf.gz; tabix -f ${data}p${i}_common.vcf.gz; done

for i in {1..3}; do gunzip -f -k ${data}p${i}_common.vcf.gz > ${data}p${i}_common.vcf ; done

bcftools merge ${data}p1_common.vcf.gz ${data}p2_common.vcf.gz  --force-samples -o ${data}p_donor.vcf

mv ${data}p3_common.vcf ${data}p_target.vcf

for i in `seq 1 40`; do echo "p1:$i p1 1" > ${data}popnames.ids; done
for i in `seq 1 40`; do echo "p2:$i p2 1" >> ${data}popnames.ids; done
for i in `seq 1 1000`; do echo "$i p_target 1" >> ${data}popnames.ids; done

# The above generates the required input to the software
# The below begins running the software, with the input data of p_donor.vcf, p_target.vcf and popnames.ids

# we need the map showing the recombination distance
$plink --vcf ${data}p_donor.vcf --recode --out ${data}p

# Longmatchquery -- find matches longer than L.
${bindir}./lmq_1swp_dpbwt$exe -i ${data}p_donor.vcf -q ${data}p_target.vcf -m -L 1 -o ${matches}target_match.txt
${bindir}./lmq_1swp_dpbwt$exe -i ${data}p_donor.vcf -q ${data}p_donor.vcf -m -L 1 -o ${matches}donor_match.txt

# split it into a match file for a target/donor individual
seq 0 999 | parallel --jobs 10 "grep -w q{} matches/target_match.txt > matches/target_match{}.txt"
seq 0 79 | parallel --jobs 10 "grep -w q{} matches/donor_match.txt > matches/donor_match{}.txt"
