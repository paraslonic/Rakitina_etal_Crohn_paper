library("gplots")

source("ecoliGroups.r")

gene.table = read.delim("ortho_table.txt", check.names = F)
rownames(gene.table) = paste(gene.table$id, gene.table$product, sep = "___")
gene.table = gene.table[,-c(1,2,match("strains", colnames(gene.table)))]
write.csv(gene.table, "genes_raw.csv")
gene.bool = gene.table
gene.bool[gene.bool > 1] = 1

strains = colnames(gene.bool)

crohns.all = union(crohns.all,strains[grep("GUTCD",strains)])
commensal = setdiff(strains, crohns.all)

# - -----------------------------------------------------------------------
# CALC   P VALUE   (FISHER) -----------------------------------------------

# get group counts
crohnall.columns = (colnames(gene.bool) %in% crohns.all)
crohnall.yes = apply(gene.bool[, crohnall.columns], 1, sum)
crohnall.no = sum(crohnall.columns) -  crohnall.yes

commensal.columns = (colnames(gene.bool) %in% commensal)
commensal.yes = apply(gene.bool[, commensal.columns], 1, sum)
commensal.no = sum(commensal.columns) - commensal.yes

# fisher
calcP = function(groups) { return(apply(groups, 1, function(x) { m1 = matrix(x,2,2); ft = fisher.test(m1); return (ft$p.value) } )) }

### commensal
bytypes.commensal = cbind(crohnall.yes, crohnall.no, commensal.yes, commensal.no)
pvalues.commensal = calcP(bytypes.commensal) 
pvalues.adj.commensal = p.adjust(pvalues.commensal)
pvalues.adj.commensal_bonferroni = p.adjust(pvalues.commensal, method = "bonferroni")

comb = data.frame(bytypes.commensal, pvalues.commensal,pvalues.adj.commensal)
write.csv(comb, "Suppltable-s03-4-5.csv", quote = FALSE)

pdf("Figure7.pdf")
the.gene.bool = gene.bool[which(pvalues.commensal < 0.01),]
the.gene.bool. = the.gene.bool

heatmap.2(as.matrix(the.gene.bool.), col = c("gray20","orange"), trace = "none",
          margins=c(8,12),cexCol = 0.5, cexRow=0.5, main = "pvalue < 0.01")

the.gene.bool = gene.bool[which(pvalues.commensal < 0.05),]
the.gene.bool. = the.gene.bool
heatmap.2(as.matrix(the.gene.bool.), col = c("gray20","orange"), trace = "none",
          margins=c(8,12),cexCol = 0.5, cexRow=0.1, main = "pvalue < 0.05")

dev.off()
