library("gplots")
library("vegan")

source("ecoliGroups.r")


pfam.table = read.delim("pfam_table", sep = "\t")
pfam.bool = pfam.table 
pfam.bool = pfam.bool[,-1]
pfam.bool[pfam.bool > 1] = 1

D = vegdist((pfam.bool), method = "bray")
strains = colnames(as.matrix(D))


### MDS
fit <- cmdscale(D,eig=TRUE, k=2)

x <- fit$points[,1]
y <- fit$points[,2]

# colors
rce.colors = "black"
litc.colors = "gray50"
path.colors = "pink"
nonpath.colors = rgb(0.12,0.56,1, 0.2)
commensal.colors = "darkolivegreen3"

count = length(x)
colors = rep("gray90", count)
colors[grep("RCE", strains)] = rce.colors
colors[strains %in% pathogenic ] = path.colors
colors[strains %in% healthy ] = nonpath.colors
colors[strains %in% crohns.ref ] = litc.colors
colors[strains %in% commensal ] = commensal.colors

pchs = rep(20, length(strains))
pchs[grep("RCE", strains)] = 19
pchs[strains %in% pathogenic ] = 18
pchs[strains %in% healthy ] = 15
pchs[strains %in% commensal ] = 16
pchs[strains %in% crohns.ref ] = 17



strains = sub("Escherichia_coli_", "", strains)
strains = sub("_uid.+", "", strains,perl = T)

#pdf("Figure5.pdf",)
png("Figure5.png", width = 2400, height = 2400, units = 'px', res = 300)


plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", cex.lab = 0.6, cex.axis = 0.6,
     main="", type="p", pch = pchs, xlim = c(-0.016, -0.008), ylim = c(-0.04, 0.028),
     cex = 1.5, col = colors)


legend("left",c("Crohn, this study","Crohn, other studies","Pathogenic","Non pathogenic","Commensal", "Unknown"), 
       col=  c(rce.colors, litc.colors, path.colors, nonpath.colors, commensal.colors,  "gray90")
       , cex = 1.2,  bty = "n",   y.intersp=1.2, pch=c(19,17,18,15,20))

dev.off()


