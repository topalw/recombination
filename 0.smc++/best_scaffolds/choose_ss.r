# read and format snp and N files
len <- read.table('ss.counts')
names(len) <- c('snps','ss')
ns <- read.table('reference_Ns.bed')
names(ns) <- c('ss','start','end')
ns$len <- ns$end - ns$start
nss <- aggregate(ns$len, by=list(ns$ss), FUN=sum)
# match
len$ns <- nss$x[match(len$ss,nss$Group.1)]
# plot
pdf('snps_Ns.pdf',height=8,width=12)
par(fig=c(0,1,.4,1),mar=c(0.5,.5,.5,.5))
barplot(len$snps[order(len$snps,decreasing = T)],
        names=len$ss[order(len$snps,decreasing = T)],
        cex.names=.5,las=2)
par(fig=c(0,1,0,.3),mar=c(0.5,.5,.5,.5),new=T)
barplot(len$ns[order(len$snps,decreasing = T)], 
        names=len$ss[order(len$snps,decreasing = T)],
        xaxt='n')
dev.off()

len2 <- len[len$snps > 100,]
write.table(len2$ss,'more_than_100_snps.table',quote=F,row.names = F,col.names=F)

len[len$snps < 100,]
