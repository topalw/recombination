  setwd('C:/Users/topalw//Desktop/recombination//LepMap/new_output/')
### ------ libraries -------###
library(wesanderson)


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
palette1 <- colorRampPalette(c('black','red','blue'))(13)
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
                      'n.scafs' =rep(0,1),
                      'scafs' = rep(0,1)
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
    # make temp dataframe for certain mask and LOD
    tmp.dat <- data.frame('lod'= rep(lod,maxLGs), #filled
                          'mask' = rep(mask,maxLGs), #filled
                          'LG'= seq(1:maxLGs), #filled
                          'n.scafs' = rep(NA,maxLGs), #to fill
                          'scafs' = rep(NA, maxLGs)) # to fill
    # loop through LGs and fill df
    for(lg in 1:maxLGs){
      tmp <- unique(pos$scaf[pos$LG==lg]) 
      tmp.dat$n.scafs[tmp.dat$LG==lg & 
                        tmp.dat$lod == lod & 
                          tmp.dat$mask==mask] <- length(tmp)  # n.scafs filled
      tmp.dat$scafs[tmp.dat$LG==lg 
                    & tmp.dat$lod == lod 
                      & tmp.dat$mask==mask ] <- paste(tmp,collapse = '+') #scafs filled
    }
    # update main df
    pos.dat <- rbind(pos.dat,tmp.dat) 
    }
}
# remove dummy 1st row
pos.dat <- pos.dat[2:nrow(pos.dat),]
#--------------PLOT POSITIONS ----------------#
for(mask in masks){
  # subset
  tmp.dat <- pos.dat[pos.dat$mask==mask ,]
  for(lod in unique(tmp.dat$lod)){
  # cumulative sum
  tmp.dat$cum[tmp.dat$lod==lod] <- cum.add(tmp.dat$scafs[tmp.dat$lod==lod])
  }
  plot(y=tmp.dat$cum,x=tmp.dat$LG,pch=16,col=add.alpha(palette1[tmp.dat$lod-3],0.7),
       xlim=c(0,30),ylim=c(0,30), ylab='cumulative # of uniq scaffolds',
       xlab='cumulative # of LGs', main=paste('cumulative #s  for mask',mask, sep=' '))
  abline(a=0,b=1)
  legend('topleft', pch=16,col=add.alpha(palette1,0.7),
         legend = unique(tmp.dat$lod[tmp.dat$mask==mask]),
         cex=0.7)
}

### ------------- PLOT COORDS ----------------### 

# set individual mask/ lod
mask = 0 
lod = 15
# get overview of what you are reading
pos.dat$scafs[pos.dat$mask==mask& pos.dat$lod==lod] # %in% gen.coord$scaf
# get data
positions <- read.delim(paste('infMask',mask,'/total_lod',lod,'_infMask',mask,'_map.txt.positions',
                              sep=''),h=F)
colnames(positions) <- c('scaf','pos','LG')
# set maxLGs
maxLGs <- min(n.lgs,length(unique(positions$LG)))
# create df with start finish of each LG
lg.coord <- data.frame('LG'= 1:maxLGs,
                       'start'=rep(0,maxLGs),
                       'scaf'=rep('',maxLGs),
                       'end'=rep(0,maxLGs))
for(lg in 1:maxLGs){
  tmp.start <- min(positions$pos[positions$LG==lg])
  tmp.scaf <- pos.dat$scafs[pos.dat$LG==lg & pos.dat$mask==mask &
                                   pos.dat$lod==lod]
  tmp.end <- max(positions$pos[positions$LG==lg])
  # update based on position of scaffold
  extra.coord <- min(gen.coord$cumsum[gen.coord$scaf == tmp.scaf])
  lg.coord$start[lg] <- tmp.start + extra.coord
  lg.coord$end[lg] <- tmp.end + extra.coord
  lg.coord$scaf[lg] <- tmp.scaf
}
  
plot(NA, ylim=c(0,maxLGs), xlim =c(0,max(gen.coord$cumsum)), ylab='LGs', 
     main = paste('LGs for mask', mask, 'and lod', lod, sep = ' '),xlab='gen. position')  
for(lg in 1:maxLGs){
  lines(y=rep(lg,2), x=c(lg.coord$start[lg], lg.coord$end[lg]))
}  
abline(v=gen.coord$cumsum)  
#text(y=27, x=gen.coord$cumsum, labels = gen.coord$scaf, cex = 0.4)
  
  
  