setwd('C:/Users/topalw/Desktop/PhD/Analyses/2.recombination/LepMap/new_output/')

lods <- 9:15
# naming scheme 
prefix <- '2nd_filtering/all_maf20_miss05_f.call'

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
tmp2 <- pos.dat[pos.dat$lod ==12 & pos.dat$theta == 0.03,]
# count how many scaffolds / LG
count <- data.frame('counts' = as.numeric(table(tmp2$LG)),
                    'new' = count.uniq(tmp2) 
)
goodLG <- row.names(count[which(count$new > 0 & count$counts < 4),])
# added the second of LG2 
goodLG <- append(goodLG,tmp2$LG[tmp2$scaf=='Super-Scaffold_2'][2])


# still need to add double scaff LGs like SS22+ & SS3+SS49 
# and any that bring a new LG to the situation

# that means that we choose LGs that add new scaffold information + have < 4 
# scaffolds inside. We have to manually add the 2nd SS2 LG in the list but this 
# is the only concistent split we have observed. 
# then we have to manual prune out markers from LGs that should not be there... 
# ex. in no mi maf20 miss05 dataset in lod12 the LG7 can be modified to only include SS22.
# and 100000100066  
# EXTRAS lod12-maf20 --> LG7 -- SS22&1....66; LG23; LG45 -- SS42 
#        lod13-maf20 --> LG4 -- SS22&1,...66; LG40; LG


