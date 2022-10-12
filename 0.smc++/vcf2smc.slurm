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

# run for every interval like so: sbatch script file.list.for.interval.x x

module purge
module load singularity 
export SINGULARITY_BINDPATH="/scratch,/dcsrsoft,/users,/work,/reference"

vcf="../1.1.vcfs/GR_tf1_DP_maf10_rec_${1}.vcf.gz"
smc_out="../2.1.smcfiles/d_${2}_${1}.smc.gz"
d=${2}
chr_mask="../chrmasks/${1}_mask.bed.gz"
chr=${1}

singularity run /dcsrsoft/singularity/containers/smcpp-56798c7d78b1.sif vcf2smc ${vcf} ${smc_out} -d ${d} ${d} -m ${chr_mask} ${chr} CH:`cat gr_samples.csv` 
