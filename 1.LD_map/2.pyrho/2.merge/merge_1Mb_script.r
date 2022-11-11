# Libs 
library(readr)
library(dplyr)
source('functions.r')
z <- c('ss13','ss42')

files <- fs::dir_ls(path=paste0('CH_windows/1Mb_windows/'),glob='*.table')
ch.py <- read_delim(files)
ch.py <- ch.py[! ch.py$ss %in% z, ] # no Z
ch.py <- thetaxN.to.hald(ch.py) # Hald cM & cM/Mb
# gr
files <- fs::dir_ls(path=paste0('GR_windows/1Mb_windows/'),glob='*.table')
gr.py <- read_delim(files)
gr.py <- gr.py[! gr.py$ss %in% z, ]
gr.py <- thetaxN.to.hald(gr.py)
# merge
py <- full_join(gr.py, ch.py, by=c('ss','start','end'))
names(py)[4:length(names(py))] <- c('GR_thetaxN', 'GR_cM', 'GR_cMMb',
                                    'CH_thetaxN', 'CH_cM', 'CH_cMMb')

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

py <- left_join(py, win, by=c('ss','start','end'))

write_delim(py,'merged_1Mb_dataset.table',quote = 'none')

