#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=36G
#SBATCH --job-name=FST
#SBATCH --account jgoudet_barn_owl
#SBATCH -e %x_%A_stderr
#SBATCH -o %x_%A_stdout

module load gcc bcftools bedtools2 r 

# run as sbatch script Super-Scaffold_x path/to/window/file

# make sure the window files are in the correct directory:

SCAFF=${1}
WIN=${2}

# make scaffold vcf
#bcftools view -r ${SCAFF} CH_GR_miss10.vcf.gz > ${SCAFF}.vcf
# make intersect 
bedtools intersect -b ${SCAFF}.vcf -a ${WIN} -wo | awk '{print $1, $2, $3, $5}' > ${WIN}_${SCAFF}.intersect

R --vanilla << EOF

library(hierfstat)
library(readr)
library(gaston,quietly=TRUE)
# global variables from bash
scaff <- "${SCAFF}"

# read files 
w <- as.data.frame(read_delim("$WIN", col_names=c('ss','start','end')))
v <- read.VCF(paste0(scaff,'.vcf'), BiAllelic=T, convert.chr=F)
inter <- read_delim(paste0("$WIN","_",scaff,'.intersect'), col_names = c('ss','start','end','pos'))

# after some problems with pop vector, code it manually and beg its in order 
pop <- rep(c('CH','GR'),c(92,23))

# make new columns and code dosage as matrix 
dos <- as.matrix(v)
w[,4] <- 0 
w[,5] <- NA
w[,6] <- NA
w[,7] <- NA
for(i in 1:nrow(w)){
  sites <- unlist(inter[inter[,2] == w[i,2],4])
  if(length(sites) > 0 ){
    # get index of snps
    indx <- which(v@snps[,4] %in% sites)
    # get f stats
    f <- fs.dosage(dos[,indx],pop)
    # write 
    w[i,4] <- length(sites)
    w[i,5] <- f[[4]][2,1]
    w[i,6] <- f[[4]][2,2]
    w[i,7] <- f[[4]][2,3]
    }
}

# write file 
options(scipen=9999)
write.table(w, paste0("${WIN}","_","${SCAFF}", ".Fst.table"), quote=F, row.names=F)

EOF 

