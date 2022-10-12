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
# both parameters are mandatory

# 0 - filter on 10% maf and get individual depth from file and get top 5 as distinguished
maf10_name=`echo ${1} | sed 's/.vcf.gz/_maf10.vcf.gz/'`

vcftools --gzvcf ${1} --maf 0.10 --recode --recode-INFO-all --stdout | bgzip -c > ${maf10_name}
vcftools --gzvcf ${maf10_name} --depth --stdout > ${maf10_name}.indDP
tail -n +2 ${maf10_name}.indDP | sort -k 3 -r | cut -f 1  | head -n 5 > DistinguishedIND.list

# name prep 
plink_out=`echo ${maf10_name} | sed 's/.vcf.gz/_rec/'`
rec_name=`echo ${maf10_name} | sed 's/.vcf.gz/_rec.vcf/'`

# 1 - recode missing and rename a bit

plink --vcf ${1} --keep-allele-order --allow-extra-chr --recode vcf-iid --out ${plink_out}
 
bgzip ${rec_name}
tabix ${rec_name}.gz

# 2 - read each line from scaffold file and do the split
while read -r line; 
do 

scaff_name=`echo ${rec_name} | sed "s/_rec/_rec_${line}/"`

bcftools view -r ${line} -O z -o "${scaff_name}.gz" "${rec_name}.gz"

tabix "${scaff_name}.gz"

done < ${2}



