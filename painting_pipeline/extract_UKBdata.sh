## This file is written by Daniel Lawson.

#!/bin/bash
## Run in slurm with:
## echo "./access.sh > runaccess.txt""
## sbatch_array.sh -f runaccess.txt

module load languages/java/sdk-1.8.0.141 binutils/2.26-GCCcore-5.4.0 apps/qctool/2.2.0

## Currently not used; this could allow looping over all chromosomes
chr=1

## unphased file
## bg="/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/raw_download$## phased file
bg="/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/raw_downloaded/haplotypes/ukb_hap_chr${chr}_v2.bgen"
## sample files
## Not sure what this one is but its in the wrong format
#sample="/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/sample-fil$## This is in the bgen sample format
sample="/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/id_mapping/data.chr1-22.sample"
## Where we want our final output
outroot="kb_hap_chr${chr}_v2"
tfile="$outroot.tmp.vcf.gz"
out="$outroot.vcf.gz"
## qctool options
#qctool -g $bg -threshhold 0.9 -snp-stats -osnp ukb_imp_chr${chr}_v2.snpstats -s $sample
## plink: only useful for the unphased bgen diles
# plink2 --bgen $bg ref-first --sample $sample --export vcf vcf-dosage=DS

qctool -g $bg -s $sample -og $tfile &> $tfile.log

## Convert from discretised (0/1) haplotype probabilities (GP) to phased GT format:
## https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
## These are P11, P12, P21 P22 "where Pij is the probability that haplotype i has allele j."
## So 1010 -> 0|0, 1001 -> 0|1, 0110-> 1|0 and 0101-> 1|1
zcat $tfile | head -n 1 | gzip - > $out
echo "##contig=<ID=$chr>" | gzip - >> $out
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" | gzip - >> $out
zcat $tfile | tail -n +3 | sed s'!1,0,1,0!0|0!g' | sed s'!0,1,0,1!1|1!g' | sed s'!0,1,1,0!1|0!g' | sed 's!1,0,0,1!0|1!g$
rm $tfile

## Check that the file structure is ok, and get frequency summaries
##module load apps/vcftools/0.1.17
##vcftools --gzvcf $out --freq --out $outroot

## Extract the header names
##zcat kb_hap_chr${chr}_v2.vcf.gz | head -n 3 | tail -n 1 | tr '\t' '\n' | tail -n +10 > kb_hap_chr${chr}_v2.indids.txt
