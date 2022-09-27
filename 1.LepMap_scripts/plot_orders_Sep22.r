f <- list.files('./orders/','*pruned')
# exclude merged - 
f1 <- f[-c(grep('20',f),grep('48',f))]
# exclude small lg
f1 <- f1[grep('55',f1,invert=T)]
f2 <- f[c(grep('20_gp',f),grep('48_gp',f))]
library(readr)
pdf('unmergedLGs-orders-pruned.pdf')
for(file in f1){
  d <- read_delim(paste0('./orders/',file))
  lg <- strsplit(file,split='_')[[1]][2]
  lm1 <- loess(d$m[!is.na(d$m)]~d$pos[!is.na(d$m)],span=.2)
  lm2 <- loess(d$f[!is.na(d$f)]~d$pos[!is.na(d$f)],span=.2)
  ord1 <- order(d$pos[!is.na(d$m)])
  ord2 <- order(d$pos[!is.na(d$f)])
  plot(lm1$fitted[ord1]~d$pos[!is.na(d$m)][ord1],type='l',
       ylab='cM',xlab='Position along the LG',
       main=paste('LG',lg,'-',unique(d$ss)),lwd=2,
       ylim=c(0,max(c(d$m,d$f),na.rm=T)))
  lines(lm2$fitted[ord2]~d$pos[!is.na(d$f)][ord2],col='magenta',lwd=2)
}
dev.off()

## code and run merged 
file <- f[grep('20_gp',f)][1]
d <- read_delim(paste0('./orders/',file))
lg <- strsplit(file,split='_')[[1]][2]


d$npos <- d$pos
lg1 <- d$ss[1]
short = TRUE
for(i in 2:nrow(d)){
  if(d$ss[i] != lg1){
    index = i 
    ends <- index -1 + 
      which(cumsum(ifelse(d$ss[index:nrow(d)] == d$ss[index],0,1)) == 1)[1]
    if(ends > index + 1){
      short = FALSE
    }
    if( short ){
      d$npos[index] <- d$npos[index-1] + 1000
      d$npos[index+1] <- d$npos[index] + 1000
    } else {
      # get starting position as other ss + 1kb
      start <- d$npos[index-1] + 1000
      d$npos[index] <- start
      # calculate steps to skip for loop
      seq1 <- seq(index+1,ends-1)
      tmp <- abs(d$pos[seq1] - d$pos[seq(index,ends-2)] )
      # make everything be the start + difference
      d$npos[seq1] <- rep(start,length(seq1)) + tmp
      # new ss first position is the last  + 1kb
      d$npos[ends] <- d$npos[ends-1] + 1000
    }
  }
}


plot(d$m[!is.na(d$m)]~d$npos[!is.na(d$m)],col=c('blue','orange')[as.factor(d$ss)],pch=16)
