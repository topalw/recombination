#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=hi
#SBATCH --account jgoudet_barn_owl

module load gcc r

R --vanilla << EOF 
	options(scipen=9999)
	library(readr)
	files <- fs::dir_ls(path='CH_m256_n92_pyrho_res/CH_windows/1kb_windows/',glob='*.table')
	d <- as.data.frame(read_delim(files))
	d[,5] <- NA
	names(d)[5] <- 'hi'
	for(scaff in unique(d[,1])){
		indexes <- which(d[,1]== scaff)
		for(i in indexes){
	        # get surrounding but not outside range and not focal
			range <- seq(i-40,i+40)
			range <- range[range %in% indexes & range != i]
			d[i,5] <- d[i,4]/mean(d[range,4],na.rm=T)
		}
	}
	write.table(d,'CH_1kb_total_hi.table',quote=F,row.names=F)
EOF


