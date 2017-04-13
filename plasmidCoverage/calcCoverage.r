library("IRanges")

calcCov <- function(path, ref.length){
  files = dir(path)
  
  cov.list =  list()
  for (f in files){
    if(file.size(paste0(path,f)) == 0) { cov.list[[f]]=0; next }
    t = read.delim(paste0(path,f), head = FALSE)
    ranges = reduce(IRanges(start = t$V7, end = t$V8))
    sum.length = sum(width(ranges))
    coverage = sum.length/ref.length
    cov.list[[f]]=coverage
  }
  
  cov = data.frame(strains = names(cov.list),cov = unlist(cov.list))
  return(cov)
}

cov_jj = calcCov("blastout_jj/", 55956)
cov_lf = calcCov("blastout_lf/", 108379)

par(mar=c(5.1,12,4.1,2.1))
pdf("plasmid_coverage.pdf")
barplot(cov_lf$cov, names.arg = cov_lf$strains, las=2, col = "dodgerblue3", 
        cex.names=  0.6, hor=TRUE, xlab="sequence coverage", main = "pLF82")

barplot(cov_jj$cov, names.arg = cov_jj$strains, las=2, col = "dodgerblue3", 
        cex.names=  0.6, hor=TRUE, xlab="sequence coverage", main = "pJJ1886_4")

dev.off()

t = data.frame(cov_lf$cov , cov_jj$cov)
rownames(t) = cov_lf$strains
colnames(t) = c("pLF82","pJJ1886_4")
write.csv(t,"plasmid_coverage.csv", row.names=FALSE)
