library("IRanges")

path = "blastout/"
ref.length = 108379
files = dir(path)

cov.list =  list()
for (f in files){
  t = read.delim(paste0(path,f), head = FALSE)
  ranges = reduce(IRanges(start = t$V7, end = t$V8))
  sum.length = sum(width(ranges))
  coverage = sum.length/ref.length
  cov.list[[f]]=coverage
}

cov = unlist(cov.list)
pdf("plf82_coverage.pdf")
par(mar=c(5.1,12,4.1,2.1))
barplot(cov, las=2, col = "dodgerblue3", cex.names=  0.6, hor=TRUE, xlab="pLF82 coverage")
dev.off()
cov = data.frame(names(cov),cov)
colnames(cov) = c("strain","plasmid coverage")
write.csv(cov,"coverage.csv", row.names=FALSE)
