
pop = args[1]

ns <- read.table('reference_Ns.bed')
colnames(ns) <- c('ss',
                  'start','end')
ns$len <- ns$end - ns$start
# anything less than 500 is not considered 
ns <- ns[ns$len >= 500,]
# Read file and rm first empty line
d <- read.table(paste0(pop,'_1000_allScaffs.table'))
d <- d[-1,]
# prune 1 chr as example
for(scaff in unique(ns$ss)){
  tmp.n <- ns[ns$ss == scaff, ]
  for(i in 1:nrow(tmp.n)){
    start <- as.numeric(tmp.n$start[i]) - 500
    windows <- which(d$start[d$scaff==scaff] %in% seq(start,tmp.n$end[i]))
    d <- d[-windows,]
  }
}

write.table(d, paste0(pop,'_1000_allScaffs_N50.table'), quote = F, row.names = F)