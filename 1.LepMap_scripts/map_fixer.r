
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Usage is Rscript map_fixer.r path/to/snps.txt path/to/.legend path/to/map.txt", 
       call.=FALSE)
} 
# read files
snps <- read.delim(args[1]) # CHR\tPOS
lgnd <- read.csv(args[2]) # lg,scaf
map <- read.delim(args[3]) # shitty header + 1 column
colnames(map) <- 'lg'
# sanity check
if(nrow(snps) != nrow(map)){
  stop('Lengths of map and snps do not match')
}
# make sure you only have the scaffs you want / lg
for(lg in lgnd$lg){
  indxs <- which(map$lg == lg)
  list.of.scaffs <- strsplit(lgnd$scaf[lgnd$lg==lg],split=';')[[1]]
  # fix positions where scaffold of snps and legend do not match for a given LG
  # AKA remove snps of weird scaffold from LGs we keep 
  map$lg[indxs[which(! snps$CHR[indxs] %in% list.of.scaffs)]] <- 0
  
}
# remove positions that do not match the LGs in the legend
map$lg[which(! map$lg %in% lgnd$lg)] <- 0 

colnames(map) <- '#lg' # make comment line 

write.table(map,paste(args[3],'_new.txt',sep=''),quote = F, row.names = F)

