#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=smcplot
#SBATCH --account jgoudet_barn_owl
#SBATCH -e TEMP/%j_stderr
#SBATCH -o TEMP/%j_stdout


# 1 is the output name  
# 2 is the input file(s)

module purge
module load singularity
export SINGULARITY_BINDPATH="/scratch,/dcsrsoft,/users,/work,/reference"


singularity run /dcsrsoft/singularity/containers/smcpp-56798c7d78b1.sif plot -c ${1} ${2}
