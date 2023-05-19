#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=30:0:0
#SBATCH --mem=10G
#SBATCH --array=1-9
#SBATCH --account=sscm011962

module load languages/gcc/10.4.0

chr=21

case ${SLURM_ARRAY_TASK_ID} in
1)
  pop='british1'
  ;;
2)
  pop='british2'
  ;;
3)
  pop='british3'
  ;;
4)
  pop='indian'
  ;;
5)
  pop='irish'
  ;;
6)
  pop='caribbean'
  ;;
7)
  pop='african'
  ;;
8)
  pop='pakistani'
  ;;
9)
  pop='chinese'
  ;;
esac

/usr/bin/time -v ./../test.exe -donorfile chr${chr}_1000G.phase.gz -targetfile chr${chr}_UKB_${pop}.phase.gz -targetname ../${pop}.txt -paintingfile painting_chr${chr}_${pop}.txt.gz -LDAfile LDA_chr${chr}_${pop}.txt.gz -LDASfile LDAS_chr${chr}_${pop}.txt.gz -mapfile chr${chr}map.txt -popfile ../popnames.txt
