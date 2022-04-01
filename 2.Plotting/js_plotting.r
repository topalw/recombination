### ------ libraries -------###
library(wesanderson)
library(ggplot2)
library(ggpubr)

###---------------FUNCTIONS-----------------###
# alpha function - <3 Tristan <3 #
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
#cumulative sum function 
cum.add <- function(scafs, separator='+'){
  # scafs vector should be 'separator' collapsed names of scaffolds
  tmp.list <- NULL
  cum.sum <- rep(0,length(scafs))
  # for 1st element
  found.scafs <- strsplit(scafs[1],split=separator, fixed = T)[[1]]
  tmp.list <- c(tmp.list,found.scafs)
  cum.sum[1] <- length(tmp.list)
  for(i in 2:length(scafs)){
    found.scafs <- strsplit(scafs[i],split=separator, fixed = T)[[1]]
    # which are new
    new <- unique(found.scafs[which(! found.scafs %in% tmp.list)])
    # append new
    tmp.list <- c(tmp.list, new)
    # update sum
    cum.sum[i] <- cum.sum[i-1] + length(new)
  }
  return(cum.sum)
}

### -------------- DATA -------------###

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (with/path/prefix) minLOD maxLOD.n", call.=FALSE)
} else if (length(args) == 1 ) {
  # default output file
  args[2] = 9
  args[3] = 15
}

lods <- seq(args[2],args[3])
print(lods)
# naming scheme 
prefix <- args[1] 
#example prefix
#'/users/atopalou/scratch/1.recombination/3_LepMap_results/all_maf20_miss05_f.call_lod12_theta0.03_map.txt_new.txt_joinSingles_'

### ------------ POSITIONS ---------------- ###

# make df to fill
pos.dat <- data.frame('lod'= rep(0,1),
                      'theta' = rep(0,1),
                      'LG' = rep(0,1),
                      'scaf' = rep(0,1),
                      'markers'=rep(0,1)
)

### FOR LOOP FOR POSITIONS ###

# lod range
# real distribution of markers 
real.dist <- read.csv(paste(strsplit(prefix,split='_lod')[[1]][1],'_dist.csv',sep=''),h=F)
colnames(real.dist) <- c('markers', 'scaf')
real.dist <- real.dist[order(real.dist$markers, decreasing = TRUE),]
real.dist$LG <- 1:nrow(real.dist)

# how many LGs to consider? 
n.lgs <-nrow(real.dist)
n.lgs <- 50
# no reason to consider more than the # of scafs in vcf?

for( theta in c(0.03)){
  # loop through lods
for(lod in lods){
  # read positions file  
  pos <- read.delim(h=F, file = 
          paste(prefix,lod,'.txt.positions'
                sep=''))
  colnames(pos) <- c('scaf','pos','LG')
  # get max LGs 
  maxLGs <- min(length(unique(pos$LG)), n.lgs)
  # loop through LGs and fill df
  for(lg in 1:maxLGs){
    scafs.in.lg <- unique(pos$scaf[pos$LG==lg])
    # more than 1 scaf / LG
      rep.num <- length(scafs.in.lg)
      markers.list <- rep(0, rep.num)
      # get number of markers / scaffold / lG
      for(i in 1:rep.num){
        markers.list[i] <- length(pos$pos[pos$LG == lg & 
                                            pos$scaf == scafs.in.lg[i]])
      }
      tmp.dat <- data.frame('lod' = rep(lod, rep.num),
                            'theta' = rep(theta, rep.num),
                            'LG' = rep(lg, rep.num),
                            'scaf' = scafs.in.lg,
                            'markers' = markers.list)
      # update main df
      pos.dat <- rbind(pos.dat,tmp.dat) 
      }
  }
}

# remove dummy 1st row
pos.dat <- pos.dat[2:nrow(pos.dat),]


#--------------PLOT GORGEOUS BARPLOT ----------------#
pos.dat$LG <- as.factor(pos.dat$LG)
#scaf.order <- order(real$markers, decreasing = T)
pos.dat$scaf <- as.factor(pos.dat$scaf)

### subset
theta = 0.03
  tmp.dat1 <- pos.dat[pos.dat$theta==theta,]
#  for(lod in unique(tmp.dat1$lod)){
  tmp.lods <- unique(tmp.dat1$lod)
  start.index <- 1
  lod = tmp.lods[start.index]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg1 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') +
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg1 <- gg1 + theme(legend.position = 'none')
    lod = tmp.lods[start.index + 1]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg2 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg2 <- gg2 + theme(legend.position = 'none')
    lod = tmp.lods[start.index + 2]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg3 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg3 <- gg3+ theme(legend.position = 'none')
#    lod = tmp.lods[start.index + 3]
#    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
#    gg4 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
#      geom_bar(stat='identity') + 
#      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
#    gg4 <- gg4 + theme(legend.position = 'none')
#    lod = tmp.lods[start.index + 4]
#    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
#    gg5 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
#      geom_bar(stat='identity') + 
#      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
#    gg5 <- gg5 + theme(legend.position = 'none')
#    lod = tmp.lods[start.index + 5]
#    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
#    gg6 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
#      geom_bar(stat='identity') + 
#      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
#    gg6 <- gg6 + theme(legend.position = 'none')
#    lod = tmp.lods[start.index + 6]
#    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
#    gg7 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
#      geom_bar(stat='identity') + 
#      theme_minimal() + 
#      ggtitle(paste('`lod=',lod, sep=''))
    gg8 <- ggplot(data=real.dist, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle('Real marker distribution')
# merge 
    fig1 <- ggarrange(gg1+ rremove('x.text')+rremove('xlab'),
              gg2+ rremove('x.text')+rremove('xlab') +rremove('ylab') ,
              gg3+ rremove('x.text')+rremove('xlab'),
             # gg4+ rremove('x.text')+rremove('xlab')+rremove('ylab'),
             # gg5+ rremove('x.text')+rremove('xlab'),
             # gg6 + rremove('x.text')+rremove('xlab')+rremove('ylab'),
             # gg7, 
	     gg8 + rremove ('ylab'),
              common.legend = TRUE, legend = 'bottom',
              ncol=2, nrow=2)

# plot 
pdf(paste(prefix,lod,'.pdf',sep=''), width=18, height=10)
title1 <- strsplit(prefix, split='/')[[1]][-1]
    
    annotate_figure(fig1, top = text_grob(paste(title1), face='bold'))

dev.off()

# FOR best lods make splitting/merging stats 
count.uniq <- function(df){
  uniq.vec <- rep('',maxLGs)
  tmp.book <- ''
  for(i in 1:maxLGs){
    lgs <- as.character(df$scaf[df$LG==i])
    new.lgs <- lgs[! lgs %in% tmp.book]
    tmp.book <- append(tmp.book, new.lgs)
    uniq.vec[i] <- length(new.lgs)
  }
  return(as.numeric(uniq.vec))
}

pal1 <- RColorBrewer::brewer.pal(length(lods),'Set2')
pdf(paste(prefix,'_diagnostics.pdf',sep=''),width=14,height=8)
for(lod in lods){
    tmp2 <- cumsum(count.uniq(tmp.dat1[tmp.dat1$lod==lod,]))
  if(lod == lods[1]){
    plot(tmp2~seq(1:maxLGs),ylab='# of new scaffolds in LG', xlab='LG', 
         main='cumulative # of Scafs in LGs', type='l',lwd=3, col=pal1[1],
         ylim=c(0,max(tmp2)))
    abline(a=0,b=1,lty=2)
  } else {
    lines(y=tmp2 ,x=seq(1:maxLGs), type='l', lwd=3, col=pal1[which(lods==lod)])
  }
}
legend('bottomright', lty=rep(1,length(lods)), lwd=rep(3,length(lods)),
       col=pal1, legend = lods, cex=0.6)
dev.off()
  
