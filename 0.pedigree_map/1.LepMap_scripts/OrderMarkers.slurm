#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --job-name=OrderMarkers
#SBATCH --account jgoudet_barn_owl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexandros.topaloudis@unil.ch

module load openjdk

### make sure there is a filename provided # 
if (( $# < 1  || $# > 3 ))
then
	echo "USAGE is sbatch OrderMarkers.slurm withPath/data withPath/map LG_to_order (1,2..N)"
	exit
fi

last=`echo ${2} | rev | cut -d '/' -f 1 | rev`
path=`echo ${2} | sed -e "s/${last}//g"`

out=${path}order_${3}.txt

### Separate Chromosomes 
java -cp ~/Lep_MAP3/bin OrderMarkers2 data=${1} map=${2} numThreads=12 chromosome=${3} usePhysical=1 > ${out}

