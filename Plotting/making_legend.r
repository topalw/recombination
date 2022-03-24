setwd('C:/Users/topalw/Desktop/PhD/Analyses/2.recombination/LepMap/new_output/')

lods <- 9:15
# naming scheme 
prefix <- '2nd_filtering/all_mi_maf10_miss05_f.call'

### ------------ FUNCTIONS ---------------- ###
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


### ------------ POSITIONS ---------------- ###

# make df to fill
pos.dat <- data.frame('lod'= rep(0,1),
                      'theta' = rep(0,1),
                      'LG' = rep(0,1),
                      'scaf' = rep(0,1),
                      'markers'=rep(0,1)
)

### FOR LOOP FOR POSITIONS ###
# real distribution of markers 
real.dist <- read.csv(paste(prefix,'_dist.csv',sep=''),h=F)
colnames(real.dist) <- c('markers', 'scaf')
real.dist <- real.dist[order(real.dist$markers, decreasing = TRUE),]
real.dist$LG <- 1:nrow(real.dist)

# how many LGs to consider? 
n.lgs <-nrow(real.dist)
n.lgs <- 200
# no reason to consider more than the # of scafs in vcf?

for( theta in c(0.03)){
  # loop through lods
for(lod in lods){
  # read positions file  
  pos <- read.delim(h=F, file = 
          paste(prefix,'_lod',lod,'_theta',theta,'_map.txt.positions',
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

# Manually curate !!!
#subset based on lod
tmp2 <- pos.dat[pos.dat$lod ==13 & pos.dat$theta == 0.03,]
# count how many scaffolds / LG
count <- data.frame('counts' = as.numeric(table(tmp2$LG)),
                    'new' = count.uniq(tmp2) 
)
# this first step will get informative LGs with a single scaffold
goodLG <- row.names(count[which(count$new > 0 & count$counts < 2),])
# added the second of LG2 
goodLG <- append(goodLG,tmp2$LG[tmp2$scaf=='Super-Scaffold_2'][2])
goodLG <- as.numeric(goodLG)[order(as.numeric(goodLG))]
goodLG
plot(cumsum(count$new)~seq(1,200))
abline(a=0,b=1)
# Others - this will be to investigate informative LGs with > 1 scaffolds
candidates <- rownames(count[which(count$new > 0 & count$counts > 1),])
tmp2[tmp2$LG %in% candidates,]
leg <- data.frame('lg' = goodLG,
                  'scaf' = tmp2$scaf[tmp2$LG %in% goodLG])

### manually add the chimeric LGs ###
# lod12-maf20 -->
# leg <- rbind(leg, data.frame('lg'=c(23,7,45),
#                              'scaf'=c('Super-Scaffold_3;Super-Scaffold_49',
#                                       'Super-Scaffold_22;Super-Scaffold_100000100066',
#                                       'Super-Scaffold_42')))
#lod13-maf20 --> 
# LG 64 doesnt have enough of SS42 to warranty forcing it 
# leg <- rbind(leg, data.frame('lg'=c(4,40),
#                              'scaf'=c('Super-Scaffold_22;Super-Scaffold_100000100066',
#                                       'Super-Scaffold_3;Super-Scaffold_49')
#                              ))
#lod12-maf10 --> 
# leg <- rbind(leg, data.frame('lg'=c(1,33),
#                              'scaf'=c('Super-Scaffold_22;Super-Scaffold_100000100066',
#                                       'Super-Scaffold_3;Super-Scaffold_49')))
#lod13-maf10 -->
# leg <- rbind(leg, data.frame('lg'=c(1),
#                              'scaf'=c('Super-Scaffold_22;Super-Scaffold_100000100066')))

write.csv(leg,paste(prefix,'_lod13_theta0.03.legend',sep=''), row.names = F, quote=F)



