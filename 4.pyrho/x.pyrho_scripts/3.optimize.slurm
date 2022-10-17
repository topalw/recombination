#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --job-name=optimize
#SBATCH --account jgoudet_barn_owl
#SBATCH -e TEMP/%x_%A_stderr
#SBATCH -o TEMP/%x_%A_stdout

# run for every SS in the vcf file

module load gcc python
source ~/rho/bin/activate

ss=${1}

vcf="../0.data/GR*${ss}.vcf.gz"
out="../3.optimize/GR_m68_${ss}"

pyrho optimize --vcffile ${vcf} --windowsize 70 --blockpenalty 15 --tablefile ../1.make_table/GR_m68 --ploidy 2 --numthreads 12 --outfile ${out}
