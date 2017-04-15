library("gplots")

source("ecoliGroups.r")

gene.table = read.delim("ortho_table.txt", check.names = F)
rownames(gene.table) = paste(gene.table$id, gene.table$product, sep = "___")
gene.table = gene.table[,-c(1,2,match("strains", colnames(gene.table)))]
gene.bool = gene.table
gene.bool[gene.bool > 1] = 1

strains = colnames(gene.bool)

crohns.all = union(crohns.all,strains[grep("GUTCD",strains)])
commensal = setdiff(strains, crohns.all)

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
pvalues = calcP(bytypes.commensal) 
pvalues.observed = pvalues

### if groups were randomly assigned?

getPvaluFromRandomized = function(){
  crohnall.columns = sample((colnames(gene.bool) %in% crohns.all))
  crohnall.yes = apply(gene.bool[, crohnall.columns], 1, sum)
  crohnall.no = sum(crohnall.columns) -  crohnall.yes
  
  commensal.columns = !crohnall.columns
  commensal.yes = apply(gene.bool[, commensal.columns], 1, sum)
  commensal.no = sum(commensal.columns) - commensal.yes
 
  calcP = function(groups) { return(apply(groups, 1, function(x) { m1 = matrix(x,2,2); ft = fisher.test(m1); return (ft$p.value) } )) }
  
  pvalues = calcP(bytypes.commensal) 
  return(pvalues)  
}

l.pval = list()
print("perfoming 1000 reshufflings...")
for(i in 1:1000){
  cat(".")
  l.pval[[i]]=getPvaluFromRandomized()
}

hist((pvalues.observed), breaks = 100, xlim=c(0,1), main = "oberved", xlab = "p-value", col = "dodgerblue")
hist((unlist(l.pval)), breaks = 50, xlim=c(0,1), main = "reshufled", xlab = "p-value", col = "dodgerblue")

pval.rand = unlist(l.pval)
sum(pval.rand < 0.05)/length(l.pval)
obs.pval = sum(pvalues.observed < 0.05)
pvalue = sum(sapply(l.pval, function(x) {sum(x < 0.05)}) >= obs.pval)/length(l.pval)
print("pvalue=",pvalue)
