library(readr)
library(dplyr)
source('functions.r')

files <- fs::dir_ls(path='forced_maps/', glob='*.perf')

orders <- data.frame('ss'=NA,'pos'=NA,'m'=NA,'f'=NA,'lg'=NA)

for(f in files){
  d <- read_delim(f, col_names = c('ss','pos','m','f'))
  lg <- as.numeric(which(files == f))
  d$ss <- rename.ss(d$ss)
  d$m <- prune.ends(d$m,0.1,2)
  d$f <- prune.ends(d$f,0.1,2)
  d <- d[! is.na(d$m) & ! is.na(d$f),]
  d$lg <- rep(lg, nrow(d))
  orders <- rbind(orders,d)
}
orders <- orders[2:nrow(orders),]

write_delim(orders,'full_forced_orders.txt',quote='none')
