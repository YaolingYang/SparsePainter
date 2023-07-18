#!/bin/bash

bgzip p0.vcf
bgzip p1.vcf
bgzip p3.vcf

bcftools index p0.vcf.gz
bcftools index p1.vcf.gz
bcftools index p3.vcf.gz

//find common SNPs and output them in output_dir/sites.txt
bcftools isec -n=3 -p output_dir p0.vcf.gz p1.vcf.gz p3.vcf.gz

bcftools view -R output_dir/sites.txt -Ov -o p0_common.vcf p0.vcf.gz
bcftools view -R output_dir/sites.txt -Ov -o p1_common.vcf p1.vcf.gz
bcftools view -R output_dir/sites.txt -Ov -o sim_target.vcf p3.vcf.gz

bgzip p0_common.vcf
bgzip p1_common.vcf

bcftools index p0_common.vcf.gz
bcftools index p1_common.vcf.gz

bcftools merge p0_common.vcf.gz p1_common.vcf.gz -Oz -o sim_ref.vcf --force-samples

//create sim_map.txt with fixed recombination rate of 5e-8 Morgan/b
echo -e "pd\tgd" > sim_map.txt
awk '{print $2"\t"5e-8*$2}' output_dir/sites.txt >> sim_map.txt

//create sim_popnames.txt
bcftools query -l sim_ref.vcf > refsamples.txt
awk '{if (NR <= 500) print $1"\t0"; else print $1"\t1"}' refsamples.txt > sim_popnames.txt

//extract target samples and create sim_targetname.txt
bcftools query -l sim_target.vcf > sim_targetname.txt

//convert vcf to phase format
pbwt -readVcfGT sim_ref.vcf -writePhase sim_ref.phase
pbwt -readVcfGT sim_target.vcf -writePhase sim_target.phase

//run HNPaint
./HMPaint.exe -reffile sim_ref.phase -targetfile sim_target.phase -popfile sim_popnames.txt -mapfile sim_map.txt -targetname sim_targetname.txt -L_minmatch 10 -matchfrac 0.01 -AAS 0 -out sim_2pop
