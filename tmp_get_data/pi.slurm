#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=pi
#SBATCH --account jgoudet_barn_owl
#SBATCH -o %x_%A.stdout
#SBTACH -e %x_%A.stderr

module load gcc samtools bcftools vcftools 

# 1 is pop 
# 2 is super-scaffold awk from list

vcf="CH_GR_miss10.vcf.gz"
pop=${1}

vcftools --gzvcf ${vcf} \
--window-pi 100000 --window-pi-step 100000 \
--chr ${2} --keep unrel_lists/${pop}_unrelated_0.03125.list \
--out ${pop}_${2}_100kb




