module load languages/r/4.1.0
module load libs/gsl/1.16-gcc-9.1.0
module load languages/gcc/9.1.0
module load apps/bcftools-1.9-74/1.9-74
module load apps/samtools/1.13-30
module load apps/vcftools/0.1.17

slim small_sim_bbpp.txt

vcftools --vcf p2.vcf --positions-overlap p1.vcf --out p2new --recode
mv p2new.recode.vcf p2new.recode.tmp.vcf
vcftools --vcf p2new.recode.tmp.vcf --positions-overlap p3.vcf --out p2new --recode

vcftools --vcf p1.vcf --positions-overlap p2new.recode.vcf --out p1new --recode
vcftools --vcf p3.vcf --positions-overlap p2new.recode.vcf --out p3new --recode

for i in `seq 1 3`; do mv p"$i"new.recode.vcf p"$i"new.vcf; done

for i in {1..3}; do bgzip -k -f p${i}new.vcf; tabix p${i}new.vcf.gz; done

for i in {1..3}; do bcftools norm -d all p${i}new.vcf.gz | bgzip -f > p${i}new.bak.vcf.gz; mv p${i}new.bak.vcf.gz p${i}new.vcf.gz; tabix -f p${i}new.vcf.gz; done

for i in {1..3}; do zcat p${i}new.vcf.gz > p${i}new.vcf ; done

bgzip -f p1new.vcf
bgzip -f p2new.vcf
bgzip -f p3new.vcf

bcftools reheader -s names_p1.txt p1new.vcf.gz -o p1_new.vcf.gz
bcftools reheader -s names_p2.txt p2new.vcf.gz -o p2_new.vcf.gz

mv p3new.vcf.gz p3.vcf.gz

tabix p1_new.vcf.gz
tabix p2_new.vcf.gz
tabix p3.vcf.gz

bcftools merge -0 p1_new.vcf.gz p2_new.vcf.gz p3.vcf.gz -o p_merged_raw.vcf
bcftools view -q 0.01:minor p_merged_raw.vcf -o p_merged_maffiltered.vcf
cat p_merged_maffiltered.vcf | sed 's!0/0!0|0!g' > p_merged_filtered.vcf

vcf2cp.pl -m 2000 -J p_merged_filtered.vcf p_merged_filtered
makeuniformrecfile.pl -c 0.000001 p_merged_filtered.phase p_merged_filtered.rec

grep popNeinf p_merged_filtered.cp

# assume N_e is 90.6536
for i in `seq 1 500`; do
fs cp -b -t popnames.ids -f popnames.donor $i $i  -r p_merged_filtered.rec -n 90.6536 -M 0 -k 1 -g p_merged_filtered.phase -o p_merged_filtered/stage7/test_p_merged_filtered_stage7_tmp_mainrun.linked_file1_ind$i; done

Rscript simulation.R
