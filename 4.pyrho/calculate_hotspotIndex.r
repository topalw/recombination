args = commandArgs(trailingOnly=TRUE)
pop = args[1]

d <- read.table(paste0('1kb_tables/',pop,'_1000_allScaffs_N50.table'),h=T)
colnames(d) <- c('start','end','theta','ss')
# find relative theta to background
d$hi <- rep(0,nrow(d))

for(scaff in unique(d$ss)){
        indexes <- which(d$ss == scaff)
        for(i in indexes){
                # get surrounding but not outside range and not focal
                range <- seq(i-40,i+40)
                range <- range[range %in% indexes & range != i]
                d$hi[i] <- d$theta[i]/mean(d$theta[range])
        }
}

write.table(d,paste0('1kb_tables/',pop,'_1000_allScaffs_N50_hi.table'),quote=F,row.names=F)
