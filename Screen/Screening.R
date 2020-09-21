library(data.table)
library(Rcmdr)
library(GOplot)
gse54129_ranked_genes <- read.csv("~/Documents/GitHub/cxxc5_rotation_project/Screen/gse54129_ranked_genes.txt", sep="")
gse54129_degs <- subset(gse54129_ranked_genes, subset=adj.P.Val < 0.05 & 
                          +                             abs(logFC) > 1.0)
gse54129_up <- subset(gse54129_degs, subset=logFC > 1.0)
gse54129_up_missing <- subset(gse54129_up, subset=Gene.symbol == "")
gse54129_up_cleaned <- subset(gse54129_up, subset=Gene.symbol != "")

gse54129_down <- subset(gse54129_degs, subset = logFC > -1.0)
gse54129_down_missing <- subset(gse54129_down, subset=Gene.symbol == "")
gse54129_down_cleaned <- subset(gse54129_down, subset=Gene.symbol != "")


gse66229_ranked_genes <- read.csv("~/Documents/GitHub/cxxc5_rotation_project/Screen/gse66229_ranked_genes.txt", sep="")
gse66229_degs <- subset(gse66229_ranked_genes, subset=adj.P.Val < 0.05 & 
                          +                             abs(logFC) > 1.0)
gse66229_up <- subset(gse66229_degs, subset=logFC > 1.0)
gse66229_up_missing <- subset(gse66229_up, subset=Gene.symbol == "")
gse66229_up_cleaned <- subset(gse66229_up, subset=Gene.symbol != "")

gse66229_down <- subset(gse66229_degs, subset = logFC < -1.0)
gse66229_down_missing <- subset(gse66229_down, subset=Gene.symbol == "")
gse66229_down_cleaned <- subset(gse66229_down, subset=Gene.symbol != "")


gse19826_ranked_genes <- read.csv("~/Documents/GitHub/cxxc5_rotation_project/Screen/gse19826_ranked_genes.txt", sep="")
gse19826_degs <- subset(gse19826_ranked_genes, subset=adj.P.Val < 0.05 & 
                          +                             abs(logFC) > 1.0)
gse19826_up <- subset(gse19826_degs, subset=logFC > 1.0)
gse19826_up_missing <- subset(gse19826_up, subset=Gene.symbol == "")
gse19826_up_cleaned <- subset(gse19826_up, subset=Gene.symbol != "")

gse19826_down <- subset(gse19826_degs, subset = logFC < -1.0)
gse19826_down_missing <- subset(gse19826_down, subset=Gene.symbol == "")
gse19826_down_cleaned <- subset(gse19826_down, subset=Gene.symbol != "")

#Need to take boxplots of each dataset (from GEO2R) showing normalization of expression values
#Need to Create Venn Diagrams of up/downregulated genes from each dataset
data1 <- subset(gse54129_degs, select=c(ID, logFC))
data2 <- subset(gse66229_degs, select=c(ID,logFC))
data3 <- subset(gse19826_degs, select=c(ID, logFC))

write.csv(data1, file = "Degs1.csv")
write.csv(data2, file = "Degs2.csv")
write.csv(data3, file = "Degs3.csv")
GOVenn(data1, data2, data3, label=c("GSE54129", "GSE66229", "GSE19826"))

quartz.save("venn_diagram.png")

overlap <- c(overlap_genes$V1)
overlapData1 <- data1[data1$ID %in% overlap,]
overlapData2 <- data2[data2$ID %in% overlap,]
overlapData3 <- data3[data3$ID %in% overlap,]

data1_com_up <- subset(overlapData1, subset=logFC >1.0)
data1_com_down <- subset(overlapData1, subset=logFC < -1.0)

data2_com_up <- subset(overlapData2, subset=logFC >1.0)
data2_com_down <- subset(overlapData2, subset=logFC < -1.0)

data3_com_up <- subset(overlapData3, subset=logFC >1.0)
data3_com_down <- subset(overlapData3, subset=logFC < -1.0)


write.csv(overlapData1, file = "Genes1.csv")
write.csv(overlapData2, file = "Genes2.csv")
write.csv(overlapData3, file = "Genes3.csv")


a <- subset(gse54129_degs, subset=Gene.symbol %in% overlap)
b <- subset(gse66229_degs, subset=Gene.symbol %in% overlap)
c <- subset(gse19826_degs, subset=Gene.symbol %in% overlap)


write.csv(a, file = "Genes1_ID.csv")
write.csv(b, file = "Genes2_ID.csv")
write.csv(c, file = "Genes3_ID.csv")


