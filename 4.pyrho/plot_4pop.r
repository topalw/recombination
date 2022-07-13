library(RColorBrewer)
pops <- c('ch','is','pt','gr')
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
col.codex <- data.frame('pop'=pops,
                        'col'= add.alpha(brewer.pal(4,'Dark2'),.6))
scaffs <- c("Super-Scaffold_10","Super-Scaffold_100000100052","Super-Scaffold_100000100064",
            "Super-Scaffold_100000100078","Super-Scaffold_100000100237","Super-Scaffold_100000100381",
            "Super-Scaffold_1000006","Super-Scaffold_12","Super-Scaffold_14","Super-Scaffold_15",
            "Super-Scaffold_16","Super-Scaffold_17","Super-Scaffold_18","Super-Scaffold_19",
            "Super-Scaffold_2","Super-Scaffold_20","Super-Scaffold_200000100656",
            "Super-Scaffold_200000105744","Super-Scaffold_20000014","Super-Scaffold_200000146",
            "Super-Scaffold_20000015","Super-Scaffold_200000175","Super-Scaffold_200000178","Super-Scaffold_2000002059","Super-Scaffold_20000042",
            "Super-Scaffold_2000008","Super-Scaffold_21","Super-Scaffold_22","Super-Scaffold_23","Super-Scaffold_24",
            "Super-Scaffold_25","Super-Scaffold_26","Super-Scaffold_27","Super-Scaffold_28","Super-Scaffold_29","Super-Scaffold_3",
            "Super-Scaffold_30","Super-Scaffold_31","Super-Scaffold_32","Super-Scaffold_33",
            "Super-Scaffold_34","Super-Scaffold_35","Super-Scaffold_36","Super-Scaffold_37",
            "Super-Scaffold_38","Super-Scaffold_39","Super-Scaffold_4","Super-Scaffold_40",
            "Super-Scaffold_41","Super-Scaffold_43","Super-Scaffold_44","Super-Scaffold_45",
            "Super-Scaffold_46","Super-Scaffold_47","Super-Scaffold_48","Super-Scaffold_49",
            "Super-Scaffold_5","Super-Scaffold_6","Super-Scaffold_7","Super-Scaffold_8","Super-Scaffold_9")
step.b <- 10000
for(scaff in scaffs){
  pdf(paste0(scaff,'_4pop_',step.b,'.pdf'))
  # sanity checkc 
  files <- list.files('.',pattern=paste0(scaff,'$'))
  if(length(files)==length(pops)){
for(pop in pops){
#reading pyrho file
dat <- read.table(paste0(pop,'_m126_mu46e9_',scaff))
colnames(dat) <- c('start','end','theta')
# make length 
dat$length <- dat$end - dat$start
dat$total <- dat$length*dat$theta
dat$mid <- (dat$end + dat$start)/2
# 1 kb windows 
a <- range(c(dat$start,dat$end))[1]
z <- range(c(dat$start,dat$end))[2]
kb <- data.frame('start'=seq(a,(z - step.b),by=step.b),'end'=seq(a+step.b,z,by=step.b))
kb$tot <- rep(0,nrow(kb))
# TO DO : MAKE FUNCTION THAT TAKES WINDOWS AND DAT 
for(i in 1:nrow(kb)){
    # first the start interval  
    first <- which(dat$start <= kb$start[i] & dat$end > kb$start[i])
    # finding the end interval
    final <- which(dat$start < kb$end[i] & dat$end >= kb$end[i])
    # check if only one 
    if(final - first == 0){
      kb$tot[i] <- dat$theta[first]*1e3
    # check if only 2 (no whole segments) so theta.first * (distance to end) +
    # theta.final * (distance to start)
    } else if (final - first == 1) {
      tot.first <- dat$theta[first] * (dat$end[first] - kb$start[i])
      tot.last <- dat$theta[last] * (kb$end[i] - dat$start[final])
      kb$tot[i] <- tot.first + tot.last
    } else {
    # add whole segments
    # segments between start and finish are all included so len*theta
    whole <- sum(dat$total[(first+1):(final-1)])
    # add the first and last segments 
    tot.first <- dat$theta[first] * (dat$end[first] - kb$start[i])
    tot.last <- dat$theta[last] * (kb$end[i] - dat$start[final])
    kb$tot[i] <- whole + tot.first + tot.last
    } } 
if(pop == pops[1]){
plot(y=kb$tot,x=seq(1:nrow(kb)),type='l',col=col.codex$col[match(pop,col.codex$pop)]
     ,ylab='sum recombination frequency',
     xlab=paste('position in ',step.b,' bp',sep=''),
     main=paste(scaff,sep=' '),ylim=c(0,6e-3), lty=which(pops==pop),lwd=2)
} else {
  lines(y=kb$tot,x=seq(1:nrow(kb)),col=col.codex$col[match(pop,col.codex$pop)],
        lty=which(pops==pop),lwd=2)
}
}
legend('topright',lty=seq(1:4),lwd=rep(2,4),col=col.codex$col,legend=col.codex$pop)
dev.off()
}
}