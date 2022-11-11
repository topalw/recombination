#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --job-name=GetData
#SBATCH --account jgoudet_barn_owl
#SBATCH -o %x_%A.stdout
#SBTACH -e %x_%A.stderr

module load gcc bcftools vcftools htslib

pop=${1}
files=`ls ~/common_work/3000owls/4_REF_PANEL/data/filteredVCFs/Renamed_TechFiltersMask_indDP*.gz`
#echo $files

bcftools concat $files | \
bcftools reheader -s newnames.list - |\
bcftools view -S unrel_lists/${pop}_unrelated_0.03125.list -v snps -m 2 -M 2 - | \
vcftools --vcf - --recode --recode-INFO-all --max-missing 0.9 --stdout | bgzip -c > ${pop}_miss10.vcf.gz 

tabix ${pop}_miss10.vcf.gz


