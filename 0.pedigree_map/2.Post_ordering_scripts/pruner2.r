library(readr)
library(dplyr)
source('functions.r')

files <- fs::dir_ls(path='forced_orders/', glob='*.perf')

orders <- data.frame('ss'=character(),'pos'=integer(),
                     'm'=double(),'f'=double(),'lg'=integer())

for(f in files){
  d <- read_delim(f, col_names = c('ss','pos','m','f'))
  lg <- as.numeric(which(files == f))
  d$ss <- rename.ss(d$ss)
  d$m <- prune.ends(d$m,0.1,2)
  d$f <- prune.ends(d$f,0.1,2)
  d <- d[! is.na(d$m) & ! is.na(d$f),]
  d$lg <- rep(lg, nrow(d))
  orders <- rbind(orders,d)
}
write_delim(orders,'full_forced_orders.txt',quote='none')


# Make individual LG marey maps and write sex-specific loess file 
p <- viridis::viridis(3) 
for(lg in unique(orders$lg)){
  pdf(paste0('forced_marey_lg',lg,'.pdf'))
  lmm <- loess(orders$m[orders$lg == lg]~orders$pos[orders$lg == lg],span=0.2)
  lmf <- loess(orders$f[orders$lg == lg]~orders$pos[orders$lg == lg],span=0.2)
  write_delim(data.frame('fitted'=lmm$fitted, 'pos'=lmm$x),paste0('lg_',lg,'loess_male'))
  write_delim(data.frame('fitted'=lmf$fitted, 'pos'=lmf$x),paste0('lg_',lg,'loess_female'))
  maxcm <- unique(max(c(orders$m[orders$lg == lg],orders$f[orders$lg == lg]),na.rm=T))
  plot(lmm$fitted[order(lmm$x)]~lmm$x[order(lmm$x)], type='l',col=p[1],lwd=2,
       xlab='position along linkage group', ylab='cM',ylim=c(0,maxcm),
       main=paste0('LG ', lg))
  lines(lmf$fitted[order(lmf$x)]~lmf$x[order(lmf$x)], col=p[2],lwd=2)
  legend('bottomright',lwd=c(2,2),col=p[c(1,2)],legend=c('Males', 'Females'))
  print(lg)
  dev.off()
}

# read all sex-specific files and merge them 
files <- fs::dir_ls('.',glob='*ale')

d <- data.frame('fitted_m'=double(),
                'pos'=integer(),
                'fitted_f'=double(),
                'lg'=integer(),
                'ss'=character())
for(lg in unique(orders$lg)){
  m <- read_delim(paste0('lg_',lg,'loess_male'))
  f <- read_delim(paste0('lg_',lg,'loess_female'))
  tmp <- dplyr::full_join(m,f,by='pos')
  names(tmp) <- c('fitted_m','pos','fitted_f')
  tmp$lg <- rep(lg,nrow(tmp))
  tmp$ss <- rep(unique(orders$ss[orders$lg == lg])[1],nrow(tmp))
  d <- rbind(d,tmp)
}
readr::write_delim(d,'forced_loess_male_female.df')

# Make sex-averaged loess dataset and write to file 
orders$sex.av <- (orders$m + orders$f) / 2
dav <- data.frame('fitted_av'=double(),
                  'pos'=integer(),
                  'lg'=integer(),
                  'ss'=character())
for(lg in unique(orders$lg)){
  lmav <- loess(orders$sex.av[orders$lg == lg] ~ orders$pos[orders$lg == lg],span=0.2)
  tmp <- data.frame('fitted_av'= lmav$fitted,
                    'pos'= lmav$x,
                    'lg'= rep(lg,length(lmav$x)),
                    'ss'= rep(unique(orders$ss[orders$lg == lg])[1],length(lmav$x)))
  dav <- rbind(dav,tmp)
  print(which(unique(orders$lg)==lg))
}
names(dav) <- c('fitted_av','pos','lg','ss')
readr::write_delim(dav,'forced_loess_sexaveraged.df')


