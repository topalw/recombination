#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --job-name=JoinSing
#SBATCH --account jgoudet_barn_owl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexandros.topaloudis@unil.ch

module load openjdk

start=`date +%s`

### make sure there is a filename provided # 
if (( $# < 1  || $# > 3 ))
then
	echo "USAGE is sbatch JoinSingles.slurm withPath/File withPath/map LODscore"
	exit
fi

lod=${3:-10}
out="${2}_joinSingles_${lod}.txt"
### Separate Chromosomes 
java -cp ~/Lep_MAP3/bin JoinSingles2All data=${1} lodLimit=${lod} map=${2} numThreads=48 > ${out}

### make counts and fix it 
fixedfile="${out}.counts.fixed" 

# since JS2 has more columns take the first one and delete after ifs
cut -f 1 ${out} > tmp_list.txt

if [ -e  $fixedfile ] 
then 
	echo "Fixed positions found"
else
	sort tmp_list.txt | uniq -c | sort > "${out}.counts" 
	tail -n +2 "${out}.counts" | sed -e 's/^ *//g' > ${fixedfile}
	rm "${out}.counts"
fi

### make positions file
snpsfile="${1}.snps.txt"
if [ -e ${snpsfile} ]
then 
	paste "${1}.snps.txt" tmp_list.txt | awk '($3>0)' > "${out}.positions"
else
	awk '(NR>=7)' "${1}"_f.call | cut -f 1,2 > ${snpsfile} 
	paste "${1}.snps.txt" tmp_list.txt | awk '($3>0)' > "${out}.positions"
fi
# delete
rm tmp_list.txt

end=`date +%s`
# this is never printing lol
echo "(($end-$start))"	
