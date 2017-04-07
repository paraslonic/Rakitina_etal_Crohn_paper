library("gplots")
library("proxy")
library("ape")
#library("MASS")

setwd("/data4/bio/runs-manolov/ecoli_crohn/analyse/pfamProfile")
bytypes = read.delim("results/bytypes_combined.csv")
head(pfam.bool)
pfam.table = read.delim("results//pfam_table")
pfam.bool = pfam.table 
pfam.bool[pfam.bool > 1] = 1
#pfam.bool = pfam.bool[-c(4,13),]

#pfam.bool = t(pfam.bool)

D = dist(pfam.bool, method = "binary")
strains = colnames(as.matrix(D))

### MDS
fit <- cmdscale(D,eig=TRUE, k=2) # k is the number of dim
#fit <- isoMDS(D, k=2) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]

#pdf("mds_bray.pdf")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="MDS. (bray distance)", type="points", col = colors, pch = 16)
legend("top", legend = c("Crohn","Crohn(ref)","non pathogenic", "others"),box.lwd=0,
       fill = c(Colors[["crohn"]],Colors[["crohn.ref"]], Colors[["healthy"]],"gray"),
       cex = 0.7)
#dev.off()



###############
#### C O L O R S   

Colors = list()
Colors[["crohn"]] = "red"
Colors[["crohn.ref"]] = "orange"
Colors[["healthy"]] = "dark green"

colors = rep("gray", length(strains))
crohns = c("RCE01", "RCE02","RCE03", "RCE04","RCE05", "RCE06","RCE07", "RCE08","RCE09","RCE10", "RCE11")
crohns.ref = c("Escherichia_coli_LF82_uid161965", "Escherichia_coli_UM146_uid162043","HM605","Escherichia_coli_O83_H1_NRG_857C_uid161987")
crohns.all = c(crohns, crohns.ref)

healthy = c("B354","Escherichia_coli__BL21_Gold_DE3_pLysS_AG__uid59245","Escherichia_coli_ATCC_8739_uid58783","Escherichia_coli_B_REL606_uid58803","Escherichia_coli_BL21_DE3__uid161947","Escherichia_coli_BL21_DE3__uid161949","Escherichia_coli_BW2952_uid59391","Escherichia_coli_DH1_uid161951","Escherichia_coli_DH1_uid162051","Escherichia_coli_ED1a_uid59379","Escherichia_coli_HS_uid58393","Escherichia_coli_IAI1_uid59377","Escherichia_coli_K_12_substr__DH10B_uid58979","Escherichia_coli_K_12_substr__MDS42_uid193705","Escherichia_coli_K_12_substr__MG1655_uid57779","Escherichia_coli_K_12_substr__W3110_uid161931","Escherichia_coli_KO11FL_uid162099","Escherichia_coli_KO11FL_uid52593","Escherichia_coli_LY180_uid219461","Escherichia_coli_SE11_uid59425","Escherichia_coli_SE15_uid161939","Escherichia_coli_W_uid162011","Escherichia_coli_W_uid162101","KLY","KTE142","KTE233","MC4100","Nissle_1917")

colors[ strains %in% crohns  ] = Colors[["crohn"]]
colors[ strains %in% crohns.ref  ] = Colors[["crohn.ref"]]
colors[ strains %in% healthy  ] = Colors[["healthy"]]




