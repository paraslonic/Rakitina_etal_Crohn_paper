library("stringr")

groups_genes = read.delim("ortho_table_names.txt",head = TRUE)
groups_genes = subset(groups_genes, groups_genes$LF82 != "")

og_id = data.frame(og = as.character(groups_genes$id), 
                   id = as.character(groups_genes$LF82), stringsAsFactors = FALSE)

### split by commas
.og_id = subset(og_id, !grepl(",", og_id$id))
for (i in grep(",", og_id$id)){
  for(x in str_split(og_id$id[i],",")[[1]]){
    .og_id <<- rbind(.og_id, c(og_id$og[i], x))
  }
}
og_id = .og_id

t = fread("table.txt")
og_pval = data.frame() 
og = sapply(t$id, function(x) { str_split(x,"___")[[1]][1] } )
pval = t$pvalues.commensal

og_pval = data.frame(og, pval)

T = merge(og_pval, og_id, by="og")
l = lapply(T$id, function(x) { str_split(x,"\\|")[[1]][4:6]} )
T = data.frame(T, do.call(rbind, l))
colnames(T) = c("og","pval","id", "contig", "start","stop")
T$start = as.numeric(as.character(T$start))
T$stop = as.numeric(as.character(T$stop))
T$pos =  0.5*(T$start + T$stop)

write.table(T, "ongenome.tab", quote = FALSE, sep="\t", row.names=FALSE)


T.chr = subset(T, T$contig == "NC_011993")
pdf("pval_on_genome.pdf")
plot(1,1,xlim=c(0,max(T.chr$pos)),ylim=c(0,1), type = "n")
for(i in 1:length(T.chr$pos)){
  p = T.chr$pval[i]
  ypos = 0.8
  if(T.chr$pval[i] <= 0.05) ypos = 0.9
  lines(c(T.chr$pos[i],T.chr$pos[i]), c(0.2,ypos),
        col = rgb(0,p,p,0.5),cex = 0.001)
}
dev.off()