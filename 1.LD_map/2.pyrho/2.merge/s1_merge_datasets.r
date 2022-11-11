## MERGER MONDAY  31.10 ##

# Libs 
library(readr)
library(dplyr)
source('functions.r')
# Global parameters  
w <- 100 # window size in number of kb
hotspot = F
lepmap = F
if(w == 1 ){ hotspot=T }
if(w == 100){ lepmap =T }
z <- c('ss13','ss42')

# 1. PYRHO 

# ch 
files <- fs::dir_ls(path=paste0('CH_windows/',w,'kb_windows/'),glob='*.table')
ch.py <- read_delim(files)
ch.py <- ch.py[! ch.py$ss %in% z, ] # no Z
ch.py <- thetaxN.to.hald(ch.py) # Hald cM & cM/Mb
# gr
files <- fs::dir_ls(path=paste0('GR_windows/',w,'kb_windows/'),glob='*.table')
gr.py <- read_delim(files)
gr.py <- gr.py[! gr.py$ss %in% z, ]
gr.py <- thetaxN.to.hald(gr.py)
# merge
py <- full_join(gr.py, ch.py, by=c('ss','start','end'))
names(py)[4:length(names(py))] <- c('GR_thetaxN', 'GR_cM', 'GR_cMMb',
                                    'CH_thetaxN', 'CH_cM', 'CH_cMMb')
# 2. Pi 

# ch
files <- fs::dir_ls(path=paste0('CH_pi/pi_',w,'kb/'),glob='*.pi')
ch.pi <- read_delim(files)[,c(1,3,4,5)]
names(ch.pi) <- c('ss','end','CH_nvar_pi','CH_pi')
ch.pi$ss <- rename.ss(ch.pi$ss) # rename Super-Scaffold_x to ssx
py <- left_join(py,ch.pi, by=c('ss','end'))
# gr
files <- fs::dir_ls(path=paste0('GR_pi/pi_',w,'kb/'),glob='*.pi')
gr.pi <- read_delim(files)[,c(1,3,4,5)]
names(gr.pi) <- c('ss','end','GR_nvar_pi','GR_pi')
gr.pi$ss <- rename.ss(gr.pi$ss)
py <- left_join(py,gr.pi, by=c('ss','end'))


# 3. Fst 

files <- fs::dir_ls(path=paste0(w,'kb_Fst'),glob='*.pruned')
fst <- read_delim(files,
                  col_names = c('ss','start','end','nvar','CH_pfst','GR_pfst','fst'))
fst$ss <- rename.ss(fst$ss)
py <- left_join(py,fst,by=c('ss','start','end'))

# 4. GC
gc <- read_delim(paste0('GC_windows/',w,'kb_GC.table'))[,c(1,2,3,5,10)]
names(gc) <- c('ss','start','end','gc','n')
gc$ss <- rename.ss(gc$ss)
py <- left_join(py,gc,by=c('ss','start','end'))

# 5. TSS
tss <- read_delim(paste0('Gene_windows/',w,'kb_tss.table'),
                  col_names = c('ss','start','end','tss'))
tss$ss <- rename.ss(tss$ss)
py <- left_join(py,tss,by=c('ss','start','end'))

# 6. TEs 
tes <- read_delim(paste0('TE_windows/',w,'kb_tes.table'),
                  col_names = c('ss','start','end','te','te_start',
                                'te_end','overlap'))
tes$ss <- rename.ss(tes$ss)
tes2 <- aggregate(tes$overlap, by=list(tes$ss,tes$start,tes$end), FUN=sum)
names(tes2) <- c('ss','start','end','te_sum')
py <- left_join(py, tes2, by=c('ss','start','end'))

# 7. introns & exons 
introns <- read_delim(paste0('Gene_windows/',w,'kb_introns.table'),
                      col_names=c('ss','start','end','ss2','ss2_start',
                                  'ss2_end','overlap'))
introns$ss <- rename.ss(introns$ss)
intr2 <- aggregate(introns$overlap, by=list(introns$ss,introns$start,introns$end),
                   FUN=sum)
names(intr2) <- c('ss','start','end','intron_sum')

py <- left_join(py, intr2, by=c('ss','start','end'))
exons <- read_delim(paste0('Gene_windows/',w,'kb_exons.table'),
                    col_names=c('ss','start','end','ss2','ss2_start',
                                'ss2_end','overlap'))
exons$ss <- rename.ss(exons$ss)
ex2 <- aggregate(exons$overlap, by=list(exons$ss,exons$start,exons$end),
                   FUN=sum)
names(ex2) <- c('ss','start','end','exon_sum')
py <- left_join(py, ex2, by=c('ss','start','end'))

# 8. conditional hotspots 
if(hotspot){
  ch.hot <- read_delim(paste0('CH_windows/CH_1kb_total_hi.table'))[,c(1,2,3,5)]
  names(ch.hot)[4] <- 'ch_hi'
  py <- left_join(py, ch.hot, by=c('ss','start','end'))
  gr.hot <- read_delim(paste0('GR_windows/GR_1kb_total_hi.table'))[,c(1,2,3,5)]
  names(gr.hot)[4] <- 'gr_hi'
  py <- left_join(py, gr.hot, by=c('ss','start','end'))
}
# 9. conditional LepMap
if(lepmap){
  lo <- read_delim('forced_loess_sexaveraged.df')
  lo.sex <- read_delim('forced_loess_male_female.df')
  win <- py[,c(1,2,3)]
  win$ped.av <- rep(NA,nrow(win))
  win$ped.m <- rep(NA,nrow(win))
  win$ped.f <- rep(NA,nrow(win))
  for(ss in unique(win$ss)){
    w1 <- win[win$ss==ss,]
    l <- lo[lo$ss == ss,]
    lsex <- lo.sex[lo.sex$ss == ss,]
    win$ped.av[win$ss == ss] <- window.ped(w1,l)
    win[win$ss == ss, c(5,6)]  <- window.ped(w1,lsex,sex.based = T)
  }
}
py <- left_join(py, win, by=c('ss','start','end'))

write_delim(py,paste0('merged_',w,'kb_dataset.table'),quote = 'none')
