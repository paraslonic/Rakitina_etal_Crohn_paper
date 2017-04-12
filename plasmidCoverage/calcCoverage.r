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
heatmap.2(t(as.matrix(t)))

par(mar=c(18,8,4,2))
image(1:nrow(t), 1:ncol(t), as.matrix(t), axes=FALSE, 
      col = colorRampPalette(c("white","blue"))(100), xlab="",ylab="")
axis(1, 1:nrow(t), rownames(t), cex.axis = 0.8, las=3)
axis(2, 1:ncol(t), c("pLF82","pJJ1886_4"), cex.axis = 0.8, las=1,col = "white", tcl = 0)
grid(nrow(t),ncol(t),lwd = 1,lty=1, col = "gray50")

colnames(t) = c("pLF82","pJJ1886_4")
heatmap.2(t(as.matrix(t)), col = colorRampPalette(c("white","blue"))(100),
          margins=c(18,12), dendrogram = "none", cexRow = 2, )

pdf("pJJ1886_4_coverage.pdf")
par(mar=c(5.1,12,4.1,2.1))
barplot(cov, las=2, col = "dodgerblue3", cex.names=  0.6, hor=TRUE, xlab="pJJ1886_4 coverage")
dev.off()

colnames(cov) = c("strain","plasmid coverage")
write.csv(cov,"coverage_pJJ1886_4.csv", row.names=FALSE)
