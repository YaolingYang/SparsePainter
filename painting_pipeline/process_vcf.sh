##painting pipelineL: paint UKB using 1000G, and investigate the difference of LDAS between populations

##run in /mnt/storage/scratch/ip21972/1000GUKB/chr${chr}

##Required files:

##1)kb_hap_chr${chr}_v2.vcf.gz  --UKBdata  
## cp ../../processUKB/kb_hap_chr${chr}_v2.vcf.gz .

##2)chr${chr}.1kg.phase3.v5a.vcf.gz ---downloaded from https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf 
##  wget https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr${chr}.1kg.phase3.v5a.vcf.gz

##3)refidx.txt  --- match file of donor samples and ancestries, downloaded from https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel

##4)chr${chr}.map --- downloaded from https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps

##5)popnames.txt
##awk '{print $1, $2}' refidx.txt | awk '!seen[$2]++ {print $2}' | sort > temp_ancestries.txt
##awk '{print $1"\t"NR-1}' temp_ancestries.txt > ancestry_to_number.txt 
##awk 'NR==FNR {a[$1]=$2; next} {print $1"\t"a[$2]}' ancestry_to_number.txt <(awk '{print $1, $2}' refidx.txt) > popnames.txt
##rm temp_ancestries.txt

##6)samples_all.txt, british1.txt, british2.txt, british3.txt, irish.txt, indian.txt, caribbean.txt, african.txt, pakistani.txt, chinese.txt

##7)test.exe

##8)beagle.jar

##9)pbwt should be downloaded and put into ~/bin

##10)getmap.R

module load languages/r/4.2.1
module load languages/gcc/10.4.0
module load apps/vcftools/0.1.17
module load apps/bcftools-1.9-74/1.9-74
module load libs/htslib/1.12

chr=16

cp ../../processUKB/kb_hap_chr${chr}_v2.vcf.gz .
wget https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr${chr}.1kg.phase3.v5a.vcf.gz

## Remove SNPs with multi alleles, multi rs names and MAF < 0.5% in both files
bcftools view -Oz -o chr${chr}_UKB_filtered_temp2.vcf.gz kb_hap_chr${chr}_v2.vcf.gz --exclude 'ID=="." || STRLEN(REF)>1 || STRLEN(ALT)>1 || N_ALT>1 || MAF<0.005'
bcftools index chr${chr}_UKB_filtered_temp2.vcf.gz
bcftools view -S ../samples_all.txt -Oz -o chr${chr}_UKB_filtered_temp.vcf.gz chr${chr}_UKB_filtered_temp2.vcf.gz
bcftools index chr${chr}.1kg.phase3.v5a.vcf.gz
bcftools view -Oz -o chr${chr}_1000G_filtered_temp.vcf.gz chr${chr}.1kg.phase3.v5a.vcf.gz --exclude 'ID=="." || STRLEN(REF)>1 || STRLEN(ALT)>1 || N_ALT>1 || MAF<0.005'
bcftools index chr${chr}_UKB_filtered_temp.vcf.gz
echo finish QC

##Get the intersection of SNPs and update vcf files
bcftools query -f'%POS\n' chr${chr}_UKB_filtered_temp.vcf.gz > UKB_positions.txt

bcftools query -f'%POS\n' chr${chr}_1000G_filtered_temp.vcf.gz > 1000G_positions.txt

comm -12 <(sort UKB_positions.txt) <(sort 1000G_positions.txt) | sort -n | awk -v chr=$chr '{print chr"\t" $1}' > shared_positions.txt

bcftools view -R shared_positions.txt chr${chr}_UKB_filtered_temp.vcf.gz -Oz -o chr${chr}_UKB_filtered.vcf.gz

bcftools index chr${chr}_1000G_filtered_temp.vcf.gz

bcftools view -R shared_positions.txt chr${chr}_1000G_filtered_temp.vcf.gz -Oz -o chr${chr}_1000G_filtered.vcf.gz

echo finish finding common SNPs

##address strand issues
bcftools index -t chr${chr}_UKB_filtered.vcf.gz

bcftools index -t chr${chr}_1000G_filtered.vcf.gz

bcftools +fixref chr${chr}_1000G_filtered.vcf.gz -Oz -o fixed_chr${chr}_1000G_filtered.vcf.gz -- -f ../hg19.fa --mode flip

bcftools +fixref chr${chr}_UKB_filtered.vcf.gz -Oz -o fixed_chr${chr}_UKB_filtered.vcf.gz -- -f ../hg19.fa --mode flip

echo finish fix ref

bcftools index -t fixed_chr${chr}_UKB_filtered.vcf.gz

bcftools index -t fixed_chr${chr}_1000G_filtered.vcf.gz

bcftools query -f '%POS\t%REF\n' fixed_chr${chr}_UKB_filtered.vcf.gz > UKB_alleles.txt
bcftools query -f '%POS\t%REF\n' fixed_chr${chr}_1000G_filtered.vcf.gz > 1000G_alleles.txt

awk -v chr=$chr 'NR==FNR {a[$1]=$2; next} ($1 in a) && ($2 == a[$1]) {print chr"\t" $1}' UKB_alleles.txt 1000G_alleles.txt > shared_positions.txt
rm chr${chr}_1000G_filtered.vcf.gz
rm chr${chr}_UKB_filtered.vcf.gz
bcftools view -R shared_positions.txt fixed_chr${chr}_1000G_filtered.vcf.gz -Oz -o chr${chr}_1000G_filtered.vcf.gz
bcftools view -R shared_positions.txt fixed_chr${chr}_UKB_filtered.vcf.gz -Oz -o chr${chr}_UKB_filtered.vcf.gz

echo finish solving strand issues

##merge these two vcf files
bcftools index -t chr${chr}_UKB_filtered.vcf.gz

bcftools index -t chr${chr}_1000G_filtered.vcf.gz

bcftools merge -Oz -o merged_chr${chr}.vcf.gz chr${chr}_UKB_filtered.vcf.gz chr${chr}_1000G_filtered.vcf.gz

echo finish merging files

##Use Beagle to phase the merged vcf file
java -Xmx100g -jar ../beagle.jar gt=merged_chr${chr}.vcf.gz map=chr${chr}.map out=merged_chr${chr}_phased

echo finish phasing

##Split the phased file into 2 files, i.e., separating donor and target samples
tabix merged_chr${chr}_phased.vcf.gz
bcftools query -l chr${chr}_1000G_filtered.vcf.gz > donorname.txt
bcftools view -S donorname.txt -Oz -o chr${chr}_1000G_phased.vcf.gz merged_chr${chr}_phased.vcf.gz
bcftools view -S ../samples_all.txt -Oz -o chr${chr}_UKB_samples.vcf.gz merged_chr${chr}_phased.vcf.gz

echo finish splitting

##Generate genetic map for SNPs
Rscript ../getmap.R ${chr}

echo finish extracting recombination map

## convert to phase files

pbwt -readVcfGT chr${chr}_1000G_phased.vcf.gz -writePhase chr${chr}_1000G.phase
gzip chr${chr}_1000G.phase

bcftools view -S ../british1.txt -Oz -o chr${chr}_UKB_british1.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_british1.vcf.gz -writePhase chr${chr}_UKB_british1.phase
gzip chr${chr}_UKB_british1.phase

bcftools view -S ../british2.txt -Oz -o chr${chr}_UKB_british2.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_british2.vcf.gz -writePhase chr${chr}_UKB_british2.phase
gzip chr${chr}_UKB_british2.phase

bcftools view -S ../british3.txt -Oz -o chr${chr}_UKB_british3.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_british3.vcf.gz -writePhase chr${chr}_UKB_british3.phase
gzip chr${chr}_UKB_british3.phase

bcftools view -S ../indian.txt -Oz -o chr${chr}_UKB_indian.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_indian.vcf.gz -writePhase chr${chr}_UKB_indian.phase
gzip chr${chr}_UKB_indian.phase

bcftools view -S ../caribbean.txt -Oz -o chr${chr}_UKB_caribbean.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_caribbean.vcf.gz -writePhase chr${chr}_UKB_caribbean.phase
gzip chr${chr}_UKB_caribbean.phase

bcftools view -S ../african.txt -Oz -o chr${chr}_UKB_african.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_african.vcf.gz -writePhase chr${chr}_UKB_african.phase
gzip chr${chr}_UKB_african.phase

bcftools view -S ../irish.txt -Oz -o chr${chr}_UKB_irish.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_irish.vcf.gz -writePhase chr${chr}_UKB_irish.phase
gzip chr${chr}_UKB_irish.phase

bcftools view -S ../pakistani.txt -Oz -o chr${chr}_UKB_pakistani.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_pakistani.vcf.gz -writePhase chr${chr}_UKB_pakistani.phase
gzip chr${chr}_UKB_pakistani.phase

bcftools view -S ../chinese.txt -Oz -o chr${chr}_UKB_chinese.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_chinese.vcf.gz -writePhase chr${chr}_UKB_chinese.phase
gzip chr${chr}_UKB_chinese.phase

echo all done

##run painting and LDA/LDAS
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_british1.phase.gz -targetname ../british1.txt -paintingfile painting_chr${chr}_british1.txt.gz -LDAfile LDA_chr${chr}_british1.txt.gz -LDASfile LDAS_chr${chr}_british1.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_british2.phase.gz -targetname ../british2.txt -paintingfile painting_chr${chr}_british2.txt.gz -LDAfile LDA_chr${chr}_british2.txt.gz -LDASfile LDAS_chr${chr}_british2.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_british3.phase.gz -targetname ../british3.txt -paintingfile painting_chr${chr}_british3.txt.gz -LDAfile LDA_chr${chr}_british3.txt.gz -LDASfile LDAS_chr${chr}_british3.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_indian.phase.gz -targetname ../indian.txt -paintingfile painting_chr${chr}_indian.txt.gz -LDAfile LDA_chr${chr}_indian.txt.gz -LDASfile LDAS_chr${chr}_indian.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_irish.phase.gz -targetname ../irish.txt -paintingfile painting_chr${chr}_irish.txt.gz -LDAfile LDA_chr${chr}_irish.txt.gz -LDASfile LDAS_chr${chr}_irish.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_caribbean.phase.gz -targetname ../caribbean.txt -paintingfile painting_chr${chr}_caribbean.txt.gz -LDAfile LDA_chr${chr}_caribbean.txt.gz -LDASfile LDAS_chr${chr}_caribbean.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_african.phase.gz -targetname ../african.txt -paintingfile painting_chr${chr}_african.txt.gz -LDAfile LDA_chr${chr}_african.txt.gz -LDASfile LDAS_chr${chr}_african.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_pakistani.phase.gz -targetname ../pakistani.txt -paintingfile painting_chr${chr}_pakistani.txt.gz -LDAfile LDA_chr${chr}_pakistani.txt.gz -LDASfile LDAS_chr${chr}_pakistani.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
##/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_chinese.phase.gz -targetname ../chinese.txt -paintingfile painting_chr${chr}_chinese.txt.gz -LDAfile LDA_chr${chr}_chinese.txt.gz -LDASfile LDAS_chr${chr}_chinese.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
