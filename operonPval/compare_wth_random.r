# Getting data

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
observed.completness  <- gene.per.op.obs$completness

# Shuffling
# because of operons in DOOR output include all genes so some operons has only one gene in it 
# we can reshuffle operon id. By doing that we preserve observed operons length distribution

l.gene.per.op.rand  <- list()

for (ind in 1:100){
  rand.op.id.v <- gene.info$OperonID
  rand.op.id.v <- rand.op.id.v[sample.int(length(rand.op.id.v))]
  gene.info$rand.op.id <- rand.op.id.v
  
  gene.per.op.rand <- merge(aggregate(data=gene.info, is.over.rep.05~rand.op.id, FUN = sum),
                       oplength,
                       by.x="rand.op.id", by.y = "OperonID")
  gene.per.op.rand$completness  <- gene.per.op.rand$is.over.rep.05/gene.per.op.rand$gene.count
  l.gene.per.op.rand[[ind]]  <- gene.per.op.rand
}

gene.per.op.rand.df = do.call(rbind, l.gene.per.op.rand)

gene.per.op.rand.max = aggregate(gene.per.op.rand.df$completness, by=list(gene.per.op.rand.df$gene.count),max)
gene.per.op.rand.mean = aggregate(gene.per.op.rand.df$completness, by=list(gene.per.op.rand.df$gene.count),mean)

gene.per.op.rand.stat = data.frame(gene.per.op.rand.mean, gene.per.op.rand.max$x)
colnames(gene.per.op.rand.stat) = c("number of genes in operon","mean fraction of genes with p-value < 0.05 (random)",
                                    "max fraction of genes with p-value < 0.05 (random)")

l.sel.operon = list()
for(i in unique(gene.per.op.rand$gene.count)){
  .obs = gene.per.op.obs[gene.per.op.obs$gene.count == i,]
 .sel =.obs$OperonID[ .obs$completness > gene.per.op.rand.max$x[gene.per.op.rand.max$Group.1 == i]]
 l.sel.operon[[i]] = .sel
}

sel.comp = subset(gene.per.op.obs, gene.per.op.obs$OperonID %in% unlist(l.sel.operon))
sel.info = subset(gene.info, gene.info$OperonID %in% unlist(l.sel.operon))

pdf("statsignificant_operons.pdf")
plot(gene.per.op.obs$gene.count, gene.per.op.obs$completness, cex.lab=1.4,
     ylab=" fraction of genes with p-value < 0.05", xlab="number of genes in operon", pch=19)

lines(gene.per.op.rand.max, col="darkolivegreen3",lty=2,lwd=2)
points(sel.comp$gene.count, sel.comp$completness, pch=19, col="dodgerblue")

dev.off()
sel.info = merge(gene.per.op.obs,sel.info,by="OperonID")
sel.info = sel.info[,c("OperonID","gene.count.x","og","prod","completness")]
colnames(sel.info)=c("OperonID","number of genes in operon","OG","product","observed fraction of genes with p-value < 0.05")
sel.info = merge(sel.info,gene.per.op.rand.stat, by="number of genes in operon")

write.table(sel.info, "statsignificant_operon_info.txt", sep="\t",row.names=F,quote=F)


