#!/usr/bin/env Rscript
# usage is Rscript 2.make_windows.r ss# 
# the intersect file is assumed to be intersect_ss#.txt 
# the new window file is assumed to be ss#_new.bed
# the output file is made from the pyrho_file
# ex ch_m126_Super-Scaffold_3 -> ch_m126_Super-Scaffold_3_w1000.table 
# for a window length of 1000

### 0 - Read arguments and sanity check 
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0 ) { 
	stop('Please specifiy the number of the super-scaffold')
}
ss <- args[1]  

### 1 - Library and file reading 
library(readr)
options(scipen=99999)
# read intersection file 
inters <- read_delim(paste0('intersect_',ss,'.txt'),
                     col_names = c('py_ss','py_start','py_end',
                                   'n_ss','n_start','n_end','overlap'))
# get window length
w <- inters$n_end[1] - inters$n_start[1] 
# read the pyrho file to get the theta from 
py <- read_delim(paste0('ch_m126_Super-Scaffold_',ss),
                 col_names = c('start','end','theta'))
# read new bed file to fill in the theta column
nf <- read_delim(paste0(ss,'_new.bed'),
                 col_names = c('ss','start','end'))

### 2 - Calculations 
# match theta to intersect file and make theta x N 
# for total prob of recombination in interval
inters$theta <- py$theta[match(inters$py_start,py$start)]
inters$thetaxN <- inters$theta*inters$overlap
# make sum of theta x N for each interval and add to nf df
tmp <- aggregate(inters$thetaxN, by=list(inters$n_start), FUN=sum)
nf$thetaxN <- tmp$x[match(nf$start,tmp$Group.1)]

### 3 - output
write_delim(nf,paste0('ch_m126_Super-Scaffold_',ss,'_w',w,'.table'),
            quote = 'none')

### END ### 
