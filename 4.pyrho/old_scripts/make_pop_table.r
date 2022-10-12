# parameters 
#setwd('C:/Users/topalw/Desktop/1.pyrho/')
pops <- c('ch','is','pt','gr')
step.b <- 1000

for(pop in pops){
  files <- list.files('.',pattern=paste0(pop,'_m126_mu46e9'))
  # make pop.df
  pop.df <- data.frame('start'=0,'end'=0,'tot'=0,'scaff'='0')
  for(file in files){
    scaff <- strsplit(file,'mu46e9_')[[1]][2]
    # sanity checkc 
      #reading pyrho file
      dat <- read.table(file)
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
          tot.last <- dat$theta[final] * (kb$end[i] - dat$start[final])
          kb$tot[i] <- tot.first + tot.last
        } else {
          # add whole segments
          # segments between start and finish are all included so len*theta
          whole <- sum(dat$total[(first+1):(final-1)])
          # add the first and last segments 
          tot.first <- dat$theta[first] * (dat$end[first] - kb$start[i])
          tot.last <- dat$theta[final] * (kb$end[i] - dat$start[final])
          kb$tot[i] <- whole + tot.first + tot.last
        } } 
      kb$scaff <- rep(scaff,nrow(kb))
      pop.df <- rbind(pop.df,kb)
  }
  write.table(pop.df,file=paste0(pop,'_',step.b,'_allScaffs.table'))

}

