# PRUNE MARKERS AT EDGE OF LGs 
library(readr)
files <- list.files('./orders','*.perf$')
# remove LGs with few loci
small.ones <- c(grep('48',files),grep('49',files),grep('52',files),grep('55',files))
files <- files[-small.ones]

# find the difference (steps in cm) in a vector
diff.all <- function(column){
  prev <- 1:(length(column)-1)
  full <- 1:length(column)
  return(column[full] - c(0,column[prev]))
}

# find max (in start) and min (in end) of steps that are > x 
find.edges <- function(column,nloc=20,cm.cut=2){
  # maximum of start
  a1 <- max(which(diff.all(column[1:nloc]) > cm.cut))
  # to escape -Inf if nothing is > x 
  a <- ifelse(is.integer(a1),a1,NA)
  # minimum of end and fix index to match length
  b <- which(diff.all(column[(length(column)-nloc):length(column)]) > cm.cut) +
    length(column) - nloc -1
  # this is due to the first element in the comparison being 0 
  # that works for the beggining but not for the end 
  # a solution would be to scale the diff to min but instead i say 
  #  if the length is > 1 (meaning it found a step other than the 1st) 
  # take the second else give NA so that the 1st element is not taken as true
  c <- ifelse(length(b) > 1, b[2],NA)
  return(c(a,c))
}

# set parameters here 
loci = 40 
cm.step = 2

for(file in files){
  d <- read_delim(paste0('orders/',file), col_types = 'cinn',
                  col_names = c('ss','pos','m','f'))
  m <- suppressWarnings(find.edges(d$m,loci,cm.step))
  f <- suppressWarnings(find.edges(d$f,loci,cm.step))
  # remove bad boys
  if(! is.na(m[1])){
    d$m[1:m[1]] <- NA
  }
  if(! is.na(m[2])){
    d$m[m[2]:length(d$m)] <- NA
  }
  if(! is.na(f[1])){
    d$f[1:f[1]] <- NA
  }
  if(! is.na(f[2])){
    d$f[f[2]:length(d$f)] <- NA
  }
  # rescale rest 
  s1 <- min(d$m,na.rm = T)
  d$m <- d$m - s1 
  s1 <- min(d$f, na.rm=T)
  d$f <- d$f - s1
  write_delim(d,paste0('orders/',file,'.',loci,'_',cm.step,'.pruned'),
              quote='none')
}


