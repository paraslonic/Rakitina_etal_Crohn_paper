# heatmap of snps on reads to assembly mapping. SNPs on heatmap are normalized to covered length, values are absolute count

library("gplots")
library("reshape2")
library("RColorBrewer")

# R E A D
L = read.delim("../readsAssembliesMap/LENGTH", head = FALSE, as.is = TRUE)
colnames(L) = c("id", "length")
L = data.table(L)
setkey(L, id)

C = read.delim("../readsAssembliesMap/COVERED", head = FALSE, as.is = TRUE)
colnames(C) = c("reads","ref", "length")
C = C[order(C$reads, C$ref),]

S = read.delim("../readsAssembliesMap/SNP_COUNT", head = FALSE, as.is = TRUE)
colnames(S) = c("reads","ref", "snp")
S = S[order(S$reads, S$ref),]

# W I D E   A N D   N O R M A L I Z E D
C.norm = C
C.norm$length = C$length/(L[C$ref]$length)

C.wide = dcast(C, reads ~ ref)
rownames(C.wide) = C.wide$reads; 
C.wide = C.wide[,-1]

C.norm.wide = dcast(C.norm, reads ~ ref)
rownames(C.norm.wide) = C.norm.wide$reads; 
C.norm.wide = C.norm.wide[,-1]

S.wide = dcast(S, reads ~ ref)
rownames(S.wide) = S.wide$reads; 
S.wide = S.wide[,-1]

S.wide.norm = S.wide/C.wide
S.wide.norm = 0.5*(S.wide.norm + t(S.wide.norm))


# P L O T
pdf("similarities.pdf")


heatmap.2(as.matrix(S.wide.norm), margins=c(8,8), trace = "none", 
          notecex=0.2,  notecol="black", cexRow = 0.7, cexCol = 0.7, 
          hclustfun = function(x) hclust(x,method = 'ward.D'),
          col = colorRampPalette(brewer.pal(10,"Blues"))(100))

dev.off()


