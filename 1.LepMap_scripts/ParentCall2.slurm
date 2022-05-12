#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --job-name=Parent_Call
#SBATCH --account jgoudet_barn_owl

module load openjdk

if (( $# < 1 || $# > 3)) 
then 
	echo "Usage is sbatch ParentCall.slurm path_to_vcf/vcf path_to_LM_ped/pedfile path_to_output/output_prefix"
	exit
fi

vcf=${1}
NFs=`echo ${vcf} | grep "/" -o | wc -l`
vcfpath1=`echo $vcf | cut -d "/" -f1-${NFs}`/

ped=${2}
out="${3:-$vcf}.call"
out2="${3:-$vcf}_f.call"

# make dataset
java -cp ~/Lep_MAP3/bin ParentCall2 halfSibs=1 data=${ped} vcfFile=${vcf} removeNonInformative=1 > ${out}
# Filter
java -cp ~/Lep_MAP3/bin Filtering2 data=${out} dataTolerance=0.01 > ${out2}

# make snps.txt
awk '(NR>=7)' ${out2} | cut -f 1,2 > ${out2}.snps.txt
tail -n +2 ${out2}.snps.txt | cut -f 1 | sort | uniq -c | sed 's/^ *//g' |sed 's/ /,/g' > ${out2}_dist.csv
