// commands to process 2-pop simulations for accuracy comparison
slim slim_sim2popsimple.txt

bgzip p0.vcf
tabix p0.vcf.gz

bgzip p1.vcf
tabix p1.vcf.gz

bgzip p3.vcf
tabix p3.vcf.gz

bcftools view -O z --exclude 'MAF<0.01' p0.vcf.gz -o p0.vcf.gz.tmp && mv p0.vcf.gz.tmp p0.vcf.gz && tabix p0.vcf.gz

bcftools view -O z --exclude 'MAF<0.01' p1.vcf.gz -o p1.vcf.gz.tmp && mv p1.vcf.gz.tmp p1.vcf.gz && tabix p1.vcf.gz

bcftools view -O z --exclude 'MAF<0.01' p3.vcf.gz -o p3.vcf.gz.tmp && mv p3.vcf.gz.tmp p3.vcf.gz && tabix p3.vcf.gz

#find common SNPs and output them in output_dir/sites.txt
bcftools isec -n=3 -p output_dir p0.vcf.gz p1.vcf.gz p3.vcf.gz

bcftools view -R output_dir/sites.txt -Ov -o p0_common.vcf p0.vcf.gz
bcftools view -R output_dir/sites.txt -Ov -o p1_common.vcf p1.vcf.gz
bcftools view -R output_dir/sites.txt -Ov -o sim_target.vcf p3.vcf.gz

#change names of sim_target.vcf
bcftools query -l sim_target.vcf > sample_names.txt
sed 's/^/p/' sample_names.txt > new_sample_names.txt
while read -r old new
do
  sed -i "s/\b$old\b/$new/g" sim_target.vcf
done < <(paste sample_names.txt new_sample_names.txt)

bgzip p0_common.vcf
bgzip p1_common.vcf

bcftools index p0_common.vcf.gz
bcftools index p1_common.vcf.gz

bcftools merge p0_common.vcf.gz p1_common.vcf.gz -Oz -o sim_ref.vcf.gz --force-samples

#create sim_map.txt with fixed recombination rate of 2e-8 Morgan/b
echo -e "pd\tgd" > sim_map.txt
awk '{print $2"\t"2e-8*$2}' output_dir/sites.txt >> sim_map.txt

#create sim_popnames.txt, change number according to the reference sizes
bcftools query -l sim_ref.vcf.gz > refsamples.txt
awk '{if (NR <= 1000) print $1"\t0"; else print $1"\t1"}' refsamples.txt > sim_popnames.txt

#extract target samples and create sim_targetname.txt
bcftools query -l sim_target.vcf > sim_targetname.txt

cp sim_ref.vcf.gz sim_ref2.vcf.gz
cp sim_target.vcf sim_target2.vcf
bgzip sim_target.vcf

#convert vcf to phase format
pbwt -readVcfGT sim_ref2.vcf.gz -writePhase sim_ref.phase
pbwt -readVcfGT sim_target2.vcf -writePhase sim_target.phase

rm sim_ref2.vcf.gz
rm sim_target2.vcf

Rscript makemap.R

java -jar flare.jar ref=sim_ref.vcf.gz ref-panel=sim_popnames.txt gt=sim_target.vcf.gz map=sim.map out=sim probs=true min-mac=1

./SparsePainter -reffile sim_ref.phase -targetfile sim_target.phase -popfile sim_popnames.txt -mapfile sim_map.txt -namefile sim_50ind.txt -L_minmatch 10 -prob -out sim_2pop

makeuniformrecfile.pl -c 2e-8 sim_target.phase recfile.rec

for i in `seq 1 50`; do 
    echo "fs cp -g sim_ref.phase -n 10 -t sim_popnames.txt -Rg sim_target.phase -r recfile.rec -f donorfile.txt $i $i -o op/sim$i -b > op/sim$i.log 2>&1" >> mycmds.txt
done
parallel -j 10 --bar < mycmds.txt
