setwd('C:/Users/topalw/Desktop/1.pyrho/')
#d <- read.delim('hyperparam/run1_extended.hyperparam')
pop <-'pt'
d <- read.table(paste0(pop,'_m126_mu46e9.hyperparam'),h=T)
d <- d[order(d$Block_Penalty,as.numeric(d$Window_Size) ),]
x <- seq(5,60,5)
y <- seq(10,100,10)
d.order <-  d[,c(1,2,3,4,5,9,10,11,12)]
colnames(d.order)
for(i in c(3:8)){
  d.order[,i] <- order(d[,match(colnames(d.order)[i],
                                colnames(d))],decreasing = T)
}  
# decreasing =F ordering for log
d.order[,9] <- order(d[,12],decreasing = F)
d.order
for(i in 1:nrow(d.order)){
  d.order$rank[i] <- sum(d.order[i,c(3:9)])
}
d.order[order(d.order$rank),c(1,2)]
top <- head(d.order[order(d.order$rank),c(1,2)],15)
top


####PLOTS####
plot(d$Spearman_Corr_100kb~d$Window_Size)
plot(d$Pearson_Corr_100kb~d$Window_Size)
# Pearson corr
z <-matrix(data=d$Pearson_Corr_1bp,byrow = T,nrow=length(x),ncol=length(y))
filled.contour(x,y,z,main='Pearson_1b',xlab='Block penalty', ylab='Window size',nlevels=40)
z <-matrix(data=d$Pearson_Corr_10kb,byrow = T,nrow=length(x),ncol=length(y))
filled.contour(x,y,z,main='Pearson_10kb',xlab='Block penalty', ylab='Window size',nlevels=40)
z <-matrix(data=d$Pearson_Corr_100kb,byrow = T,nrow=length(x),ncol=length(y))
filled.contour(x,y,z,main='Pearson_100kb',xlab='Block penalty', ylab='Window size',nlevels=40)
# L2 
z <-matrix(data=d$Log_L2,byrow = T,nrow=length(x),ncol=length(y))
filled.contour(x,y,z,main='log L2',xlab='Block penalty', ylab='Window size',nlevels=40)
# Spearman corr
z <-matrix(data=d$Spearman_Corr_1bp,byrow = T,nrow=length(x),ncol=length(y))
filled.contour(x,y,z,main='Spearman_1b',xlab='Block penalty', ylab='Window size',nlevels=40)
z <-matrix(data=d$Spearman_Corr_10kb,byrow = T,nrow=length(x),ncol=length(y))
filled.contour(x,y,z,main='Spearman_10kb',xlab='Block penalty', ylab='Window size',nlevels=40)
z <-matrix(data=d$Spearman_Corr_100kb,byrow = T,nrow=length(x),ncol=length(y))
filled.contour(x,y,z,main='Spearman_100kb',xlab='Block penalty', ylab='Window size',nlevels=40)
