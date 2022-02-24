setwd('C:/Users/topalw/Desktop/PhD/Analyses/2.recombination/LepMap/new_output/')
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
  found.scafs <- strsplit(scafs[1],split='+', fixed = T)[[1]]
  tmp.list <- c(tmp.list,found.scafs)
  cum.sum[1] <- length(tmp.list)
  for(i in 2:length(scafs)){
    found.scafs <- strsplit(scafs[i],split='+', fixed = T)[[1]]
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
# get real counts 
real <- read.delim('real.counts',sep=' ',h=F)
colnames(real) <- c('markers', 'LG')
real <- real[order(real$markers, decreasing = T),]
real$prop <- real$markers/sum(real$markers)
barplot(real$markers[ordered(real$LG)]~real$LG, xaxt='n',ylab='# of markers / scaf',
        xlab='scaffolds')
# get chromosome lengths (from markers!)
lengths <- read.delim('positions',h=F)
colnames(lengths) <- c('scaf','pos')
gen.coord <- data.frame('scaf'=NULL,'pos'=NULL)
for(scaf in unique(lengths$scaf)){
  end1 <- max(lengths$pos[lengths$scaf==scaf])
  tmp.gen.coord <- data.frame('scaf'=rep(scaf,2),
                              'pos'=c(1,end1))
  gen.coord <- rbind(gen.coord, tmp.gen.coord)
}
gen.coord$cumsum <- cumsum(gen.coord$pos)
# masks 
masks <- c(0,1,2,3)

# make palette with # of cols == # of lods
pal1.fun <- colorRampPalette(c('black','orange','blue'))
palette1 <- colorRampPalette(c('black','orange','blue'))(30)
### --------------- COUNTS --------------- ###

### FOR LOOP FOR COUNTS ###

for(mask in masks){
# read LOD values used for certain mask
lods <- read.delim(paste('infMask',mask,'/lods', sep=''),h=F)[,1]
# how many LGs should we plot?
n.lgs <- length(real$LG) + 1 # as many scaffolds+ 1
# empty df
counts.dat <- data.frame(
  'markers'=rep(0,1),
  'LG'=rep(0,1),
  'lod'=rep(0,1),
  'prop'=rep(0,1)
)
for( lod in lods ){

# ------------------- COUNTS --------------------- #    
  # read file and name columns
  dat <- read.delim(h=F, file=
    paste('infMask',mask,'/total_lod',lod,'_infMask',mask,'_map.txt.counts.fixed',
          sep=''),sep=' ')
  colnames(dat) <- c('markers','LG') 
  tots <- sum(dat[,1])
  # order by # of markers
  dat <- dat[order(dat$markers, decreasing = T),] 
  # this is going to be X (whichever of 2 is smaller)
  maxLGs <- min(nrow(dat), n.lgs)
  # choose only X rows
  dat <- dat[1:maxLGs,]
  # populate dataframe
  dat$lod <- rep(lod,nrow(dat))
  dat$prop <- dat$markers/sum(dat$markers[dat$LG!=0])
  
  counts.dat <- rbind(counts.dat, dat)
}
counts.dat <- counts.dat[2:nrow(counts.dat),]


# ------------------PLOT COUNTS----------------------- #

dat.plt <- counts.dat[counts.dat$LG!=0,]
par(mfrow=c(1,2))
plot(dat.plt$markers~dat.plt$LG, cex=1.2, 
     col = add.alpha(palette1[dat.plt$lod-3], 0.5), pch =16,
     ylab='# of markers',xlab='LG', log='y',
     main=paste('# of markers distribution / LOD for mask ',mask,sep=''))
legend('topright',pch=16,col=add.alpha(palette1,0.5),legend = lods, title='LOD', cex=.7)  

plot(dat.plt$prop~dat.plt$LG, cex=1.2, 
     col = add.alpha(palette1[dat.plt$lod-3],0.5), pch =16,
     ylab='proportion of markers',xlab='LG', log='y',
     main=paste('prop of markers for mask ',mask,sep=''),ylim=c(5e-6,8e-01))
legend('topright',pch=16,col=add.alpha(palette1,0.5),legend = lods, title='LOD', cex=.7) 
points(real$prop~seq(1,length(real$LG)), col='black',cex=1.2)
par(mfrow=c(1,1))
barplot(counts.dat$markers[counts.dat$LG==0]~counts.dat$lod[counts.dat$LG==0],
        ylab='# of singles', xlab='LOD',col=add.alpha(palette1,0.5),
        main=paste('# of markers not assigned to LGs for mask',mask, sep=''))
}

### ------------ POSITIONS ---------------- ###

# make df to fill
pos.dat <- data.frame('lod'= rep(0,1),
                      'mask' = rep(0,1),
                      'LG' = rep(0,1),
                      'scaf' = rep(0,1),
                      'markers'=rep(0,1)
)

### FOR LOOP FOR POSITIONS ###

for(mask in masks){
  # read lods file
  lods <- read.delim(paste('infMask',mask,'/lods', sep=''),h=F)[,1]
  # loop through lods
  for(lod in lods){
  # read positions file  
    pos <- read.delim(h=F, file = 
            paste('infMask',mask,'/total_lod',lod,'_infMask',mask,'_map.txt.positions',
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
                              'mask' = rep(mask, rep.num),
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
scaf.order <- order(real$markers, decreasing = T)
pos.dat$scaf <- as.factor(pos.dat$scaf)
#for(mask in masks){
  # subset
mask = 2
  tmp.dat1 <- pos.dat[pos.dat$mask==mask,]
#  for(lod in unique(tmp.dat1$lod)){
  tmp.lods <- unique(tmp.dat1$lod)
  tmp.lods
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
    lod = tmp.lods[start.index + 3]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg4 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg4 <- gg4 + theme(legend.position = 'none')
    lod = tmp.lods[start.index + 4]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg5 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg5 <- gg5 + theme(legend.position = 'none')
    lod = tmp.lods[start.index + 5]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg6 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg6 <- gg6 + theme(legend.position = 'none')
    lod = tmp.lods[start.index + 6]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg7 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg7 <- gg7 + theme(legend.position = 'none')
    lod = tmp.lods[ start.index + 7]
    tmp.dat <- tmp.dat1[tmp.dat1$lod==lod,]
    gg8 <- ggplot(data=tmp.dat, aes(x=LG, y=markers, fill=scaf)) +
      geom_bar(stat='identity') + 
      theme_minimal() + ggtitle(paste('lod=',lod, sep=''))
    gg8 <- gg8 + theme(legend.position = 'none')
    pdf(paste('mask',mask,'_gorgeous.pdf',sep=''))
    ggarrange(gg1+ rremove('x.text')+rremove('xlab'),
              gg2+ rremove('x.text')+rremove('xlab') +rremove('ylab') ,
              gg3+ rremove('x.text')+rremove('xlab'),
              gg4+ rremove('x.text')+rremove('xlab')+rremove('ylab'),
              gg5+ rremove('x.text')+rremove('xlab'),
              gg6 + rremove('x.text')+rremove('xlab')+rremove('ylab'),
              gg7 , gg8+rremove('ylab'), 
              common.legend = TRUE, legend = 'bottom',
              ncol=2, nrow=4)
    dev.off()

#}



    
    
    
    
  