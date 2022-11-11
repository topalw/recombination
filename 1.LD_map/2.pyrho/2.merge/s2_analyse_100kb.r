library(readr)
library(dplyr)
source('functions.r')
d <- read_delim('merged_100kb_dataset.table')

### 1.compare with pedigree 

# calculate total order length
orders <- read_delim('full_orders.txt')
orders$ss <- rename.ss(orders$ss)
msum <- sum(aggregate(orders$m, by=list(orders$lg), FUN=max)[,2])
fsum <- sum(aggregate(orders$f, by=list(orders$lg), FUN=max)[,2])
lepmap.length <- (fsum + msum) / 2

# corresponding length from pyrho 
sum.gr <- sum(aggregate(d$GR_cM[d$ss %in% orders$ss], 
                    by = list(d$ss[d$ss %in% orders$ss]), 
                    FUN=sum, na.rm=T)[,2])
sum.ch <- sum(aggregate(d$CH_cM[d$ss %in% orders$ss], 
                    by = list(d$ss[d$ss %in% orders$ss]), 
                    FUN=sum, na.rm=T)[,2])
# ratios
gr.ratio <-  lepmap.length / sum.gr 
ch.ratio <- lepmap.length /  sum.ch 
# scale 
d$GR_cM <- d$GR_cM * gr.ratio 
d$CH_cM <- d$CH_cM * ch.ratio

# Per Linkage group stats 

lgs <- d %>% 
      group_by(ss) %>%
      summarize(len = max(end)/1e6,
            GR_cM = sum(GR_cM, na.rm=T),
            CH_cM = sum(CH_cM, na.rm=T), 
            GR_pi = median(GR_pi, na.rm=T),
            CH_pi = median(CH_pi, na.rm=T),
            fst = median(fst, na.rm=T),
            gc = median(gc, na.rm=T),
            tss = sum(tss, na.rm =T),
            ped = sum(ped, na.rm=T),
            CH_cMMb_len = sum(CH_cM, na.rm=T)/max(end)*1e6, 
            GR_cMMb_len = sum(GR_cM, na.rm=T)/max(end)*1e6,
            ped_cMMb_len = sum(ped, na.rm=T)/max(end)*1e6)

# subset relevant ss
tmp <- lgs[lgs$ped != 0,]
# pi ~ length for microchr effect
plot(tmp$len,tmp$CH_pi,pch=16,ylab='median pi in CH',xlab='length in bp')
tmp[tmp$len==min(tmp$len),]
text(tmp$len,tmp$CH_pi,labels=tmp$ss,adj=-.5)
# cM/Mb for 
seq1 <- seq(1,100)
seq2 <-50/seq1
par(mfrow=c(1,3))
plot(tmp$len,tmp$CH_cMMb_len, pch=16, ylab = 'cM/Mb in each ss',xlab='length in Mbp')
lines(seq1,seq2)
plot(tmp$len,tmp$GR_cMMb_len, pch=16, ylab = 'cM/Mb in each ss',xlab='length in Mbp')
lines(seq1,seq2)
plot(tmp$len,tmp$ped_cMMb_len, pch=16, ylab = 'cM/Mb in each ss',xlab='length in Mbp')
lines(seq1,seq2)
dev.off()
# pi with recombination per SS
plot(tmp$CH_pi~tmp$CH_cMMb_len,pch=16, xlab= 'cM/Mb in each ss', ylab='median pi in CH')
# 
plot(tmp$CH_cMMb_len,tmp$GR_cMMb_len,pch=16,ylab='GR',xlab='CH')
abline(0,1)
plot(tmp$CH_cMMb_len,tmp$ped_cMMb_len,pch=16,ylab='ped',xlab='CH')
abline(0,1)

