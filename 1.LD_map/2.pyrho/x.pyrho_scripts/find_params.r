library(readr)
library(dplyr)
filename <- 'ch_m126.hyperparam'
d <- read_delim(filename)
# get sum of correlations
sums <- d %>% 
  select(! contains('Log') & contains('Corr'))  %>% 
  rowSums()
# add parameters
d2 <- d %>% 
  select(Block_Penalty, Window_Size)
d2$sum <- sums 
# add l2
d2$L2 <- d$L2
# plot
pdf(paste0(filename,'.pdf'))
plot(d2$L2,d2$sum,pch=16,xlab='L2',ylab='Sum of Corr')
# highlight min l2 
minl2 <- which(d2$L2==min(d2$L2))
points(d2$L2[minl2],y=d2$sum[minl2],col='red')
text(x=d2$L2[minl2],y=d2$sum[minl2],labels=paste(d2$Block_Penalty[minl2],
                               d2$Window_Size[minl2]))
# highlight max sum  
maxsum <- which(d2$sum==max(d2$sum))
points(d2$L2[maxsum],y=d2$sum[maxsum],col='red')
text(x=d2$L2[maxsum],y=d2$sum[maxsum],labels=paste(d2$Block_Penalty[maxsum],
                                                 d2$Window_Size[maxsum]))
dev.off()
