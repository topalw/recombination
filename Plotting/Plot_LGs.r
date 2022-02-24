setwd('c:/Users/topalw/Desktop/recombination/')
library(wesanderson)
positions <- read.delim('positions', h=F) # the physical positions
colnames(positions) <- c('Scaffold', 'pos') 
positions$Scaffold <- as.factor(positions$Scaffold)
truth <- table(positions$Scaffold)
truth <- truth[order(truth, decreasing = T)]
### plot size of LGs and get stats ###

# loop code should include these:
stats.df <- data.frame(matrix(data=NA,nrow=1,ncol=3))
colnames(stats.df) <-c('LGs','SNP_in_LGs','Singles')
lods <- c(3,5,10,15,20)
for(rep in lods){
  numbers <- read.delim(paste('logs/lod', rep, '.numbers',sep=''),h=F) # stats
  colnames(numbers) <- c('LGs','SNP_in_LGs','Singles')
  counts <- read.delim(paste('total_lod', rep, '.counts', sep =''), h=F, sep = ' ') 
  stats.df <- rbind(stats.df, numbers)
  #counts of snps / LG
  colnames(counts) <- c('LG','SNPs')
  counts <- counts[1:30,]
  if(rep == 3){  
    plot(counts$SNPs~counts$LG, lwd=2, col = wes_palette('IsleofDogs1')[1], type='l')
    points(y=counts$SNPs,x=counts$LG,  pch=16, cex=1.2, col = wes_palette('IsleofDogs1')[1])
  }else{
    lines(y=counts$SNPs,x=counts$LG, lwd=2,pch=16, col=wes_palette('IsleofDogs1')[which(lods == rep)])
    points(y=counts$SNPs,x=counts$LG,  pch=16, cex=1.2, col = wes_palette('IsleofDogs1')[which(lods == rep)])

  }
  abline(v=1)
}
lines(y=truth, x=1:length(truth),  pch = 16, lwd=2,col='red')
points(y=truth, x=1:length(truth),  pch = 16, cex=1.5,col='red')

legend('topright',legend=lods,fill=wes_palette('IsleofDogs1')[1:5], title = 'LOD')
# finish statistics dataframe
rownames(stats.df[2:6,]) <- lods
stats.df <- stats.df[2:6,]
stats.df$lod <- lods
print(stats.df)

### ENDS here ###
 

### Compare LGs with physical scaffolds ###