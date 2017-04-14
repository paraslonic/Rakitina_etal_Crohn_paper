system('for f in ../pfams/*.pfam; do sed -i "s/\\s\\+/\t/g" $f; done')

library("reshape2")

pfams = data.frame()
pfams.long = data.frame()

pfam.files = dir("../pfams","*.pfam")

for(f in pfam.files){ # read pfam file for each strain
  name = sub(".pfam","", f)
  t = read.delim(paste0("../pfams/",f), skip = 28, head = F)
  pfams.long = rbind(pfams.long, data.frame(t, name = paste0(t$V6, "___", t$V7, "___", t$V8)) )
  pfams = rbind(pfams, data.frame(rep(name, nrow(t)), name = paste0(t$V6, "___", t$V7, "___", t$V8))) 
}

rm(pfam.files)

colnames(pfams)=c("strain","pfam")
pfams$pfam = as.character(pfams$pfam)
pfams$strain = as.character(pfams$strain)
pfam.table = dcast(pfams, "strain~pfam", c())
rownames(pfam.table) = pfam.table$strain
pfam.table = pfam.table[,-match("strain", colnames(pfam.table))]

write.table(pfam.table, "pfam_table",quote=F, sep="\t",row.names = T)
#write.table(pfams.long, "pfams_long",quote=F, sep="\t",row.names = F)
