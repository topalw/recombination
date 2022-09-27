files <- list.files(path='./orders/',pattern='20_2.pruned')
library(readr)

orders <- data.frame('ss'=NA,'pos'=NA,'m'=NA,'f'=NA,'lg'=NA)
for(f in files){
  d <- read_delim(paste0('./orders/',f))
  lg <- strsplit(f,split='_')[[1]][2]
  d <- d[!is.na(d$m) & !is.na(d$f),]
  d$lg <- rep(lg,nrow(d))
  orders <- rbind(orders,d)
}
orders <- orders[2:nrow(orders),]
write_delim(orders,'full_orders.txt',quote='none')
