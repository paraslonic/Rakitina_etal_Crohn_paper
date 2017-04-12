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

t = data.frame(cov_lf$cov , cov_jj$cov)
rownames(t) = cov_lf$strains
colnames(t) = c("pLF82","pJJ1886_4")

pdf("pJJ1886_4_coverage.pdf")

heatmap.2(t(as.matrix(t)), col = colorRampPalette(c("white","blue"))(100),
          margins=c(18,12), dendrogram = "none", cexRow = 2)

dev.off()

write.csv(t,"plasmid_coverage.csv", row.names=FALSE)
