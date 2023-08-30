#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=30:0:0
#SBATCH --mem=10G
#SBATCH --array=1-7
#SBATCH --account=sscm011962

module load languages/gcc/10.4.0

chr=21

case ${SLURM_ARRAY_TASK_ID} in
1)
  pop='british'
  ;;
2)
  pop='indian'
  ;;
3)
  pop='irish'
  ;;
4)
  pop='caribbean'
  ;;
5)
  pop='african'
  ;;
6)
  pop='pakistani'
  ;;
7)
  pop='chinese'
  ;;
esac

/usr/bin/time -v ./../HMPaint.exe -aveSNPpaintingfile aveSNPpainting_chr${chr}_${pop}.txt.gz -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_${pop}.phase.gz -targetname ../${pop}.txt -aveindpaintingfile aveindpainting_chr${chr}_${pop}.txt.gz -LDAfile LDA_chr${chr}_${pop}.txt.gz -LDASfile LDAS_chr${chr}_${pop}.txt -AASfile AAS_chr${chr}_${pop}.txt -mapfile chr${chr}map.txt -popfile ../popnames.txt -fixlambda 164.2
