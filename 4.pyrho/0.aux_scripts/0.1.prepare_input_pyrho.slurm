#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=0.1.pyrho
#SBATCH --account jgoudet_barn_owl
#SBATCH -o %x_%A.stdout
#SBTACH -e %x_%A.stderr

module load gcc htslib vcftools 

# USAGE is sbatch 0.1.prepare_input_pyrho.slurm path_to_vcf scaffold_file 

in=`basename ${1}`
#path=`dirname ${1}`

while read -r line;
do 

out=`echo ${in} | sed "s/.vcf.gz/_mac2_thin10_${line}.vcf.gz/"`

vcftools --gzvcf ${1} --recode --recode-INFO-all --thin 10 --mac 2 --chr ${line} --stdout | bgzip > "../0.data/${out}"

tabix "../0.data/${out}"

done < ${2}


