#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=vcf2smc
#SBATCH --account jgoudet_barn_owl
#SBATCH -e TEMP/%j_stderr
#SBATCH -o TEMP/%j_stdout


module purge
module load singularity 
export SINGULARITY_BINDPATH="/scratch,/dcsrsoft,/users,/work,/reference"


# parameters needed are:
# 1 - path/to/files/XX where XX is the pop code (ex CH)
# 2-  scaffold name (ex Super-Scaffold_3 ) 
# distinguished individuals file must be named XX.distinguished where XX is the pop code and be in this directory 

pop=`basename ${1}`
ss=${2}
vcf=`ls ${1}*${ss}.vcf.gz`
chr_mask="../chrmasks/${2}_mask.bed.gz"

while read -r line; 

do

d=${line}

smc_out="../2.smcfiles/d_${d}_${pop}_${ss}"
#echo $smc_out
singularity run /dcsrsoft/singularity/containers/smcpp-56798c7d78b1.sif vcf2smc ${vcf} ${smc_out} -d ${d} ${d} -m ${chr_mask} ${ss} ${pop}:`cat ${pop}_samples.csv` 

done < ${pop}.distinguished
