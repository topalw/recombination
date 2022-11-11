args = commandArgs(trailingOnly=TRUE)
pop = args[1]

ns <- read.table('~/work/ref_genome_2020/Ns/reference_Ns.bed')
colnames(ns) <- c('ss',
                  'start','end')
ns$len <- ns$end - ns$start
# anything less than 500 is not considered
ns <- ns[ns$len >= 500,]
# Read file and rm first empty line
d <- read.table(paste0('1kb_tables/',pop,'_1000_allScaffs.table'))
d <- d[-1,]
# add scaffs that have no Ns
d1 <- d[! d$scaff %in% unique(ns$ss), ]
# prune 1 chr as example
for(scaff in unique(ns$ss)){
  tmp.n <- ns[ns$ss == scaff, ]
  tmp.d <- d[d$scaff == scaff,]
  windows <- c()
  for(i in 1:nrow(tmp.n)){
  # set range of N window where kb window would have >50% N
    start <- as.numeric(tmp.n$start[i]) - 500
    end <- tmp.n$end[i] - 499
    positions <- seq(start,end,1)
    # find those windows that start there
    windows <- append(windows, which(tmp.d$start %in% positions) )
  }
  if( length(windows) != 0 ){
  d1 <- rbind(d1, tmp.d[-windows,])
  } else { d1 <- rbind(d1, tmp.d) } # the if statement should fix the issue with scaffolds missing if no overlapping Nwindows and pyrho windows!
}
# write table
write.table(d1, paste0('1kb_tables/',pop,'_1000_allScaffs_N50.table'),quote=F,row.names=F)
