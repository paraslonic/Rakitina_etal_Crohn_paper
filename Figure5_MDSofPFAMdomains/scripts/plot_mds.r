library("gplots")
library("vegan")

source("ecoliGroups.r")


pfam.table = read.delim("pfam_table", sep = " ")
pfam.bool = pfam.table 
pfam.bool = pfam.bool[,-1]
pfam.bool[pfam.bool > 1] = 1

D = vegdist(t(pfam.bool), method = "bray")
strains = colnames(as.matrix(D))


### MDS
fit <- cmdscale(D,eig=TRUE, k=2)

x <- fit$points[,1]
y <- fit$points[,2]

# colors
rce.col = "black"
litc.col = "gray40"
path.col = "pink"
nonpath.col = "darkolivegreen3"


col = rep("gray80", length(strains))
col[grep("RCE", strains)] = rce.col
col[strains %in% pathogenic ] = path.col
col[strains %in% healthy ] = nonpath.col
col[strains %in% crohns.ref ] = litc.col

pchs = rep(20, length(strains))
pchs[grep("RCE", strains)] = 19
pchs[strains %in% pathogenic ] = 18
pchs[strains %in% healthy ] = 15
pchs[strains %in% crohns.ref ] = 17



strains = sub("Escherichia_coli_", "", strains)
strains = sub("_uid.+", "", strains,perl = T)

#pdf("Figure5.pdf",)
png("Figure5.png", width = 2400, height = 2400, units = 'px', res = 300)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", cex.lab = 0.6, cex.axis = 0.6,
     main="", type="p", col = col, pch = pchs, xlim = c(-0.041, 0.05), ylim = c(-0.05, 0.04),
     cex = 1.5)


legend("bottomright", c("Crohn, this study","Crohn, other studies","Pathogenic","Non pathogenic","Other"), col=  c(rce.col, litc.col, path.col, nonpath.col, "gray80"), cex = 1.2,  bty = "n",   y.intersp=1.2, pch=c(19,17,18,15,20))

dev.off()


