ch <- read.table('ch_1000_allScaffs_N50_hi.table',h=T)
pt <- read.table('pt_1000_allScaffs_N50_hi.table',h=T)
gr <- read.table('gr_1000_allScaffs_N50_hi.table',h=T)
is <- read.table('is_1000_allScaffs_N50_hi.table',h=T)
# make kb guide or Mb guide 
ls <- read.table('contig.lengths')
names(ls) <- c('ss','len')
# try first 100kb scale (that means only use ss that have > 50000 len)
w <- 100000 # window
big.boys <- ls$ss[ls$len > w/2]
# scaffold loop 
df.w <- data.frame('ss'=NULL,'start'=NULL,'end'=NULL,
                   'tot.ch'=NULL,'tot.pt'=NULL,
                   'tot.gr'=NULL, 'tot.is'=NULL)
for(ss in big.boys){
# make data frame
nw <- seq(0,ls$len[ls$ss==ss],by=w)
if(length(nw) %% 2 == 0){
dw <- matrix(seq(0,ls$len[ls$ss==ss],by=w),ncol=2,byrow = T)
} else { 
  dw <- matrix(c(nw,nw[length(nw)]+w),ncol=2,byrow = T)}
dw <- data.frame(dw)
names(dw) <- c('start','end')
# make columns for each pop
dw$tot.ch <- rep(NA,nrow(dw))  
dw$tot.pt <- rep(NA,nrow(dw))       
dw$tot.gr <- rep(NA,nrow(dw))       
dw$tot.is <- rep(NA,nrow(dw))     

# CH
for(i in 1:nrow(dw)){
  ch.in.w <- ch$theta[ch$ss == ss & ch$start >= dw$start[i] & ch$end < dw$end[i]]
  pt.in.w <- pt$theta[pt$ss == ss & pt$start >= dw$start[i] & pt$end < dw$end[i]]
  gr.in.w <- gr$theta[gr$ss == ss & gr$start >= dw$start[i] & gr$end < dw$end[i]]
  is.in.w <- is$theta[is$ss == ss & is$start >= dw$start[i] & is$end < dw$end[i]]
  if(length(ch.in.w) > (w/1000)/2) { 
    dw$tot.ch[i] <- sum(ch.in.w)
  } 
  if(length(pt.in.w) > (w/1000)/2) { 
    dw$tot.pt[i] <- sum(pt.in.w)
  } 
  if(length(gr.in.w) > (w/1000)/2) { 
    dw$tot.gr[i] <- sum(gr.in.w)
  } 
  if(length(is.in.w) > (w/1000)/2) { 
    dw$tot.is[i] <- sum(is.in.w)
  } 
}
# append scaffold to  df.w 
dw$ss <- rep(ss,nrow(dw))
df.w <- rbind(df.w, dw)
print(ss)
}
# REMOVE NAs and ch outliers
#d <- df.w[!is.na(df.w$tot.ch),]
#d <- d[!is.na(d$tot.pt),]
#d <- d[!is.na(d$tot.gr),]
#d <- d[!is.na(d$tot.is),]

#write.table(d,'10kb_windows_allpops.table',quote=F,row.names = F)

d <- df.w[df.w$tot.ch < 0.05,]
d <- d[! is.na(d$start),]
# order LGs based on # of windows 
ord1 <- names(table(d$ss))[order(as.numeric(table(d$ss)),decreasing = T)]
d$ss <- factor(d$ss, levels=ord1)
d <- d[order(d$ss),]
# plotting 
d$pos <- cumsum(d$start)
pops <- c('ch','pt','gr','is') 
p1 <- c('#CCB391','#8F1800','#BD6D30','#CCD9D6')
for(pop in pops){
#pdf(paste0(pop,'_100kb_map.pdf'),width=24,height=24)
par(fg='white',bg='black',col.axis='white',col.lab='white')
plot(1, type="n", xlab="Position (Mb)", ylab="theta", xlim=range(d$pos,na.rm=T),
     ylim=c(0,0.02),axes=F)
start = 1
for(ss in levels(d$ss)){
  lines(x=d$pos[d$ss==ss],
        y=d[d$ss==ss,which(pops==pop)+ 2],
        col=c(p1[match(pop,pops)],'white')[which(unique(d$ss) == ss) %% 2 +1 ],lwd=2)
}
axis(2,at=seq(0,0.02,0.005), labels=seq(0,0.02,0.005),pos=-1e9)
axis(1,at=c(seq(0,90e9,10e9),range(d$pos)[2]),labels=c(seq(0,900,by=100),round(range(d$pos)[2]/1e8)), pos=-0.0005)
#dev.off()
}
# try scaling all based on ch inferred results 
pdf('all_pops_100kb_map2.pdf',width=20,height=12)
for(pop in pops){
  ymax <-c(1,.77,.52,.27)[which(pops==pop)]
  ymin <- c(.75,.5,.25,0)[which(pops==pop)]
  margin <- c(0,0,0,2)[which(pops==pop)]
  # scale results to swiss max theta
  scaling <- sum(d$tot.ch,na.rm = T)/sum(d[,which(pops==pop)+ 2],na.rm = T)
  
  par(fg='white',bg='black',col.axis='white',col.lab='white',fig=c(0,1,ymin,ymax),new=TRUE,
      mar = c(margin,0,0,0))
  plot(1, type="n", xlab='',ylab='', xlim=range(d$pos,na.rm=T),ylim=c(0,0.015),axes=F)
  start = 1
  for(ss in levels(d$ss)){
    lines(x=d$pos[d$ss==ss],
          y=d[d$ss==ss,which(pops==pop)+ 2] * scaling,
          col=c(p1[match(pop,pops)],'darkgrey')[which(unique(d$ss) == ss) %% 2 +1 ],lwd=2)
  }
  axis(2,at=seq(0,0.015,0.005), labels=c(seq(0,0.010,0.005),''),pos=-1e9)
}
axis(1,at=c(seq(0,90e9,10e9),range(d$pos)[2]),labels=c(seq(0,900,by=100),round(range(d$pos)[2]/1e8)), pos=-0.0005)
dev.off()


 