##painting pipelineL: paint UKB using 1000G, and investigate the difference of LDAS between populations

##Required files:

##1)kb_hap_chr${chr}_v2.tmp.vcf.gz  --UKBdata

##2)chr${chr}.1kg.phase3.v5a.vcf.gz ---downloaded from https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf 
##  wget https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr21.1kg.phase3.v5a.vcf.gz

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

chr=21

## Remove SNPs with multi alleles, multi rs names and MAF < 0.5% in both files
bcftools view -Oz -o chr${chr}_UKB_filtered_temp.vcf.gz --exclude 'ID!="." && STRLEN(REF)>1 || N_ALT>1 || MAF<0.005' kb_hap_chr${chr}_v2.tmp.vcf.gz
bcftools index chr${chr}.1kg.phase3.v5a.vcf.gz
bcftools view -Oz -o chr${chr}_1000G_filtered_temp.vcf.gz --exclude 'ID!="." && STRLEN(REF)>1 || N_ALT>1 || MAF<0.005' chr${chr}.1kg.phase3.v5a.vcf.gz

##Get the intersection of SNPs and update vcf files
bcftools query -f'%POS\n' chr${chr}_UKB_filtered_temp.vcf.gz > UKB_positions.txt

bcftools query -f'%POS\n' chr${chr}_1000G_filtered_temp.vcf.gz > 1000G_positions.txt

comm -12 <(sort UKB_positions.txt) <(sort 1000G_positions.txt) | sort -n | awk -v chr="${chr}" '{print chr "\t" $1}' > shared_positions.txt

bcftools index chr${chr}_UKB_filtered_temp.vcf.gz

bcftools view -R shared_positions.txt chr${chr}_UKB_filtered_temp.vcf.gz -Oz -o chr${chr}_UKB_filtered.vcf.gz

bcftools index chr${chr}_1000G_filtered_temp.vcf.gz

bcftools view -R shared_positions.txt chr${chr}_1000G_filtered_temp.vcf.gz -Oz -o chr${chr}_1000G_filtered.vcf.gz

##merge these two vcf files
bcftools index -t chr${chr}_UKB_filtered.vcf.gz

bcftools index -t chr${chr}_1000G_filtered.vcf.gz

bcftools merge -Oz -o merged_chr${chr}.vcf.gz chr${chr}_UKB_filtered.vcf.gz chr${chr}_1000G_filtered.vcf.gz

##Use Beagle to phase the merged vcf file
java -Xmx100g -jar ../beagle.jar gt=merged_chr${chr}.vcf.gz map=chr${chr}.map out=merged_chr${chr}_phased

##Split the phased file into 2 files, i.e., separating donor and target samples
tabix merged_chr${chr}_phased.vcf.gz
bcftools query -l chr${chr}_1000G_filtered.vcf.gz > donorname.txt
bcftools view -S donorname.txt -Oz -o chr${chr}_1000G_phased.vcf.gz merged_chr${chr}_phased.vcf.gz
bcftools view -S ../samples_all.txt -Oz -o chr${chr}_UKB_samples.vcf.gz chr${chr}_UKB_phased.vcf.gz


##Generate genetic map for SNPs
Rscript ../getmap.R ${chr}

## convert to phase files

pbwt -readVcfGT chr${chr}_1000G_phased.vcf.gz -writePhase chr${chr}_1000G.phase
gzip chr${chr}_1000G.phase

bcftools view -S ../british1.txt -Oz -o chr${chr}_UKB_british1.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_british1.vcf.gz -writePhase chr${chr}_UKB_brtish1.phase
gzip chr${chr}_UKB_british1.phase

bcftools view -S ../british2.txt -Oz -o chr${chr}_UKB_british2.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_british2.vcf.gz -writePhase chr${chr}_UKB_brtish2.phase
gzip chr${chr}_UKB_british2.phase

bcftools view -S ../british3.txt -Oz -o chr${chr}_UKB_british3.vcf.gz chr${chr}_UKB_samples.vcf.gz
pbwt -readVcfGT chr${chr}_UKB_british3.vcf.gz -writePhase chr${chr}_UKB_brtish3.phase
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

##run painting and LDA/LDAS
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_british1.phase.gz -targetname british1.txt -paintingfile painting_chr${chr}_british1.txt.gz -LDAfile LDA_chr${chr}_british1.txt.gz -LDASfile LDAS_chr${chr}_british1.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_british2.phase.gz -targetname british2.txt -paintingfile painting_chr${chr}_british2.txt.gz -LDAfile LDA_chr${chr}_british2.txt.gz -LDASfile LDAS_chr${chr}_british2.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_british3.phase.gz -targetname british3.txt -paintingfile painting_chr${chr}_british3.txt.gz -LDAfile LDA_chr${chr}_british3.txt.gz -LDASfile LDAS_chr${chr}_british3.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_indian.phase.gz -targetname indian.txt -paintingfile painting_chr${chr}_indian.txt.gz -LDAfile LDA_chr${chr}_indian.txt.gz -LDASfile LDAS_chr${chr}_indian.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_irish.phase.gz -targetname irish.txt -paintingfile painting_chr${chr}_irish.txt.gz -LDAfile LDA_chr${chr}_irish.txt.gz -LDASfile LDAS_chr${chr}_irish.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_caribbean.phase.gz -targetname caribbean.txt -paintingfile painting_chr${chr}_caribbean.txt.gz -LDAfile LDA_chr${chr}_caribbean.txt.gz -LDASfile LDAS_chr${chr}_caribbean.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_african.phase.gz -targetname african.txt -paintingfile painting_chr${chr}_african.txt.gz -LDAfile LDA_chr${chr}_african.txt.gz -LDASfile LDAS_chr${chr}_african.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_pakistani.phase.gz -targetname pakistani.txt -paintingfile painting_chr${chr}_pakistani.txt.gz -LDAfile LDA_chr${chr}_pakistani.txt.gz -LDASfile LDAS_chr${chr}_pakistani.txt.gz -mapfile chr${chr}_map.txt
##/usr/bin/time -v ./test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_chinese.phase.gz -targetname chinese.txt -paintingfile painting_chr${chr}_chinese.txt.gz -LDAfile LDA_chr${chr}_chinese.txt.gz -LDASfile LDAS_chr${chr}_chinese.txt.gz -mapfile chr${chr}_map.txt
