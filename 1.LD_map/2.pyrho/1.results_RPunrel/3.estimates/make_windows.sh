#!/bin/bash

module load gcc bedtools2 r 

FILE=${1}
WIND=${2:-1000}

SS=`echo $FILE | cut -d '_' -f 4`

R --vanilla << EOF 
	options(scipen=9999)
	w <- $WIND
	library(readr)

	d <- read_delim("$FILE",
        	        col_names=c('start','end','theta'))
	ss <- strsplit("$FILE",split='Scaffold_')[[1]][2]
	scaff <- paste0('ss',ss)
	bed_py <- data.frame('ss'=rep(scaff,nrow(d)), 
	                     'start'=d[,1],
			     'end'=d[,2])
	start <- 0
	end <- (max(bed_py[,3]) %/% w ) * w
	seq1 <- seq(start,end,w) 
	seq2 <- c(seq(w, max(bed_py[,3]), w), max(bed_py[,3]))
	bed_new <- data.frame('ss' = rep(scaff,length(seq1)),
			      'start' = seq1,
			      'end' = seq2)
	write.table(bed_py, paste0(ss,'_py.bed'),quote=F,row.names=F, col.names=F,sep='\t')
	write.table(bed_new,paste0(ss,'_new.bed'),quote=F,row.names=F, col.names=F,sep='\t')
EOF

bedtools intersect -b ${SS}_new.bed -a ${SS}_py.bed -wo > intersect_${SS}.txt 

Rscript 2.make_windows.r ${SS}

rm ${SS}_new.bed ${SS}_py.bed intersect_${SS}.txt 

