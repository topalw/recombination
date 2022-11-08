### FUNCTION 1 - α 
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

### FUNCTION 2 - rename Super-Scaffold_## to ss##
rename.ss <- function(column){
  return(paste0('ss',unlist(strsplit(column,'_'))[seq(2,length(column)*2,2)]))
}

### FUNCTION 3 - make cM and cM/Mb columns from thetaxN
thetaxN.to.hald <- function(d){
  if(! 'thetaxN' %in% names(d)){
    stop('column not found')
  } else{
    # this is cM / window
    d$cM <- -.5 * log(1-2*d$thetaxN) * 100
    # depending on window we need to transform to cM/Mb
    d$cMMb <- d$cM * (1e6 /(d$end[1] - d$start[1]))
  }
  return(d)
}

### FUNCTION 4 - make windowed from orders 
open.windows <- function(windows, orders){
  width <- windows[1,3] - windows[1,2]
  windows$ped <- NA
  # verify that relevant columns exist 
  orders$av.cm <- (orders$m + orders$f)/2
  orders <- orders[! is.na(orders$av.cm), ]
  sss <- unique(windows$ss)
  sss <- sss[sss %in% unique(orders$ss)]
  for(ss in sss){
    lm1 <- loess(orders$av.cm[orders$ss==ss]~orders$pos[orders$ss==ss],span=0.2)
    for(i in 1:nrow(windows[windows$ss == ss,])){
      indx <- which(lm1$x >= windows$start[i] & lm1$x < windows$end[i])
      if(length(indx) > 0 ){
        windows$ped[windows$ss==ss][i] <- max(lm1$fitted[indx]) - min(lm1$fitted[indx]) 
      }
    }
  }
  
  return(windows)
}

### FUNCTION 5 -  remove outlier markers from LepMap3 orders 
prune.ends <- function(cMcol, perc.loc=0.05, cm.cutoff=2){
  nloc <- round(length(cMcol) * perc.loc ,0)
  # maximum of start
  a1 <- max(which(diff(cMcol[1:nloc]) > cm.cutoff)) +1 # +1 because diff prints from 2
  # to escape -Inf if nothing is > x 
  a <- ifelse(! is.na(a1 %% 1),a1,NA)
  # minimum of end and fix index to match length
  end.indx <- (length(cMcol)-nloc):length(cMcol)
  b1 <- which(diff(cMcol[end.indx]) > cm.cutoff) + length(cMcol) - nloc +1 
  # take the smallest of b1 (since it shows the first jump
  b <- ifelse(! is.na(b1 %% 1), b1, NA)[1] # escape -Inf
  a <- ifelse(! is.na(a), a, 1) # if its not NA give a if it is give the start
  b <- ifelse(! is.na(b), b, length(cMcol)) # else give the end
  # rescale everything relevat 
  cMcol[a:b] <- cMcol[a:b]  - min(cMcol[a:b])
  # set all else to NA
  true.indx <- 1:length(cMcol)
  outside <- true.indx[! true.indx %in% seq(a,b)]
  cMcol[outside] <- NA
  return(cMcol)
}