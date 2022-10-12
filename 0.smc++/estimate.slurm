#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --job-name=smcpp-est
#SBATCH --account jgoudet_barn_owl
#SBATCH -e TEMP/%j_stderr
#SBATCH -o TEMP/%j_stdout


module purge
module load singularity
export SINGULARITY_BINDPATH="/scratch,/dcsrsoft,/users,/work,/reference"
files=`ls ../2.smcfiles/*.gz`
out=${1}

singularity run /dcsrsoft/singularity/containers/smcpp-56798c7d78b1.sif estimate --spline piecewise -o ${out} --cores 12 4.6e-9 ${files} 
