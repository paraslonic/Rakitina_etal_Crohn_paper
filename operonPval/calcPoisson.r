pval <- 0.05
MIN_OPERON_LENGTH  <- 1

gene.tab  <- read.delim("lf82genes.txt", head = T)
gene.tab  <- gene.tab[,c("id","prok","prod")]

ongenome <- read.table("ongenome.tab", header = T,stringsAsFactors = F,fill = TRUE, sep = "\t", comment.char = '', quote = '')
ongenome$id <- gsub('.*\\|(.*?)\\|.*', '\\1', ongenome$id)
ongenome  <- data.frame( ongenome[,c("id","og","pval")])

gene.info  <- merge(gene.tab, ongenome)
gene.info$is.over.rep.05 <- as.numeric(gene.info$pval<pval)

opinfo  <- read.delim("door_output.txt")
opinfo <- opinfo[,c("Operon","Gene")]
colnames(opinfo) <- c("OperonID","prok")

oplength <- aggregate(opinfo$OperonID,by=list(opinfo$OperonID), FUN=length)
colnames(oplength) <- c("OperonID","gene.count")

opinfo <- merge(opinfo, oplength)

opinfo$prok <- as.character(opinfo$prok)
opinfo$prok  <- gsub(" ","",opinfo$prok)
gene.info$prok <- as.character(gene.info$prok)

gene.info <- merge(gene.info, opinfo,by="prok")

# Calc observed operon completness
gene.per.op.obs <- merge(aggregate(data=gene.info, is.over.rep.05~OperonID, FUN = sum),
                         oplength, by = "OperonID")
gene.per.op.obs$completness  <- gene.per.op.obs$is.over.rep.05/gene.per.op.obs$gene.count

### Poisson

rate = sum(ongenome$pval<0.05)/nrow(ongenome)
  
gene.per.op.obs$expected = gene.per.op.obs$gene.count*rate

pvals = apply(gene.per.op.obs,1, function(x){ppois(x[2], x[5], lower = FALSE)})
gene.per.op.obs$pvalsadj = p.adjust(pvals)
gene.per.op.obs = gene.per.op.obs[order(gene.per.op.obs$pvalsadj),]
gene.per.op.obs = gene.per.op.obs[order(gene.per.op.obs$completness, decreasing = T),]
gene.per.op.obs. = subset(gene.per.op.obs, gene.per.op.obs$pvalsadj < 0.05)
nrow(gene.per.op.obs.)
result = merge(gene.per.op.obs., gene.info, by = "OperonID")

write.table(result,"statsignificant_operon_info_Poisson.txt", quote = F, row.name = F,sep="\t")
