#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=76G
#SBATCH --job-name=maketable
#SBATCH --account jgoudet_barn_owl
#SBATCH -e TEMP/%x_%A_stderr
#SBATCH -o TEMP/%x_%A_stdout

#dcsrsoft use old 
module load gcc python 
source ~/rho/bin/activate

pyrho make_table --samplesize 186 --approx --moran_pop_size 279 \
 --numthreads 12 --mu 4.6e-9 --outfile ../1.make_table/ch_m279 \
 --smcpp_file CH.csv

