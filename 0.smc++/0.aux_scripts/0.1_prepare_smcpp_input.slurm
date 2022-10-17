#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=0.1.smcpp
#SBATCH --account jgoudet_barn_owl

module load gcc gatk vcftools bcftools samtools plink-ng/1.9.20200712
module load openjdk/1.8.0_265-b01 htslib

# run like sbatch 0.1_prepare_smcpp_input.slurm path_to_input_vcf.gz scaffold_file  
# both first and second parameters are mandatory but third defaults to CH

scaffold_file=$2
# 0 - filter on 10% missingness and get individual depth from file and get top 5 as distinguished
miss10_name=${1}
vcftools --gzvcf ${miss10_name} --depth --stdout > ${miss10_name}.indDP
tail -n +2 ${miss10_name}.indDP | sort -k 3 -r | cut -f 1  | head -n 5 > "${miss10_name}.distinguished"
# make sample csv 
tmp=`tail -n +2 "${miss10_name}.distinguished" | cut -f 1`
echo $tmp | sed 's/ /,/' > "${miss10_name}_samples.csv"

# name prep 
plink_out=`echo ${miss10_name} | sed 's/.vcf.gz/_rec/'`
rec_name=`echo ${miss10_name} | sed 's/.vcf.gz/_rec.vcf/'`

# 1 - recode missing and rename a bit

plink --vcf ${miss10_name} --keep-allele-order --allow-extra-chr --recode vcf-iid --out ${plink_out}
 
bgzip ${rec_name}
tabix ${rec_name}.gz

# 2 - read each line from scaffold file and do the split
while read -r line; 
do 

scaff_name=`echo ${rec_name} | sed "s/_rec/_rec_${line}/"`

bcftools view -r ${line} -O z -o "${scaff_name}.gz" "${rec_name}.gz"

tabix "${scaff_name}.gz"

done < ${2}



