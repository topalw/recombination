#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --job-name=SepChr1
#SBATCH --account jgoudet_barn_owl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexandros.topaloudis@unil.ch

module load openjdk

start=`date +%s`

### make sure there is a filename provided # 
if (( $# < 1  || $# > 3 ))
then
	echo "USAGE is sbatch SeparateChromosomes.slurm withPath/File LODscore mask theta"
	exit
fi

lod=${2:-10}
theta=${4:-0.03}
mask=${3:-1}
out="${1}_lod${lod}_theta${theta}_mask${mask}_map.txt"

### Separate Chromosomes 
java -cp ~/Lep_MAP3/bin SeparateChromosomes2 data=${1} lodLimit=${lod} numThreads=36 sizeLimit=3 theta=${theta} informativeMask=${mask} > ${out}

### make counts and fix it 
fixedfile="${out}.counts.fixed" 
if [ -e  $fixedfile ] 
then 
	echo "Fixed positions found"
else
	sort $out | uniq -c | sort > "${out}.counts" 
	tail -n +2 "${out}.counts" | sed -e 's/^ *//g' > ${fixedfile}
	rm "${out}.counts"
fi

### make positions file
snpsfile="${1}.snps.txt"
if [ -e ${snpsfile} ]
then 
	paste "${1}.snps.txt" "${out}" | awk '($3>0)' > "${out}.positions"
else
	awk '(NR>=7)' "${1}"_f.call | cut -f 1,2 > ${snpsfile} 
	paste "${1}.snps.txt" "${out}" | awk '($3>0)' > "${out}.positions"
fi

end=`date +%s`

echo "$((end-start))"	
