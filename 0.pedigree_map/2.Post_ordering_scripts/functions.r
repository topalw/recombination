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
window.ped <- function(windows, ped, sex.based = FALSE){
  if(sex.based){
    res <- data.frame('m'=rep(NA,nrow(windows)),
                      'f'=rep(NA,nrow(windows)))
    
  } else { res <- rep(NA,nrow(windows)) }
  for(i in 1:nrow(windows)){
    start <- windows$start[i]
    end <- windows$end[i]
    indx <- which(ped$pos < end & 
                    ped$pos >= start)
    if(length(indx) > 1 ){
    if(sex.based){
      res$m[i] <- max(ped$fitted_m[indx], na.rm=T) - 
        min(ped$fitted_m[indx], na.rm=T)
      res$f[i] <- max(ped$fitted_f[indx], na.rm=T) - 
        min(ped$fitted_f[indx], na.rm=T)
    } else {
      res[i] <- max(ped$fitted_av[indx], na.rm=T) - 
        min(ped$fitted_av[indx], na.rm=T)
    }
    } else { next }
  }
  return(res)
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