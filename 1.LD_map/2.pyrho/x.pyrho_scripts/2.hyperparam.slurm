#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --job-name=hyperparam
#SBATCH --account jgoudet_barn_owl
#SBATCH -e TEMP/%x_%A_stderr
#SBATCH -o TEMP/%x_%A_stdout

# run for every interval like so: sbatch script file.list.for.interval.x x

module load gcc python
source ~/rho/bin/activate

pyrho hyperparam --samplesize 84 --tablefile ../1.make_table/ch_m126 --mu 4.6e-9 --ploidy 2 --smcpp_file ch.csv --outfile ../2.hyperparam/ch_m126.hyperparam --numthreads 24 --windowsize 30,40,50,60,70,80,90 --blockpenalty 15,20,25,30,35,40,45,50 
