library(DGCA, quietly = TRUE)
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2, quietly= TRUE)
library(formattable)
library(data.table)


ranked_genes <- fread("gse66229_ranked_genes.txt")
top_250 <- head(ranked_genes, 250)



n_normal_samples <- 100
n_gastric_samples <- 300

cell_type <- c(rep("Normal", n_normal_samples), rep("Gastric", n_gastric_samples))
design_mat <- model.matrix(~ cell_type + 0)
colnames(design_mat) <- c("Normal", "Gastric")
x <- str(design_mat)

design_mat <- makeDesign(cell_type)
design_mat


cx1 <- "224516_s_at"
cx2 <- "233955_x_at"
cx3 <- "222996_s_at"
cxxc5_probes <- c(cx1, cx2, cx3)

top_250_genes_id <- c(top_250$ID, cxxc5_probes)

top_250_genes_exp <- gse_normal_then_gastric[gse_normal_then_gastric$X1 %in% top_250_genes_id,]
gene_names_as_rownames <- c(top_250_genes_exp$X1)
rownames(top_250_genes_exp) <- gene_names_as_rownames
top_250_genes_exp <- subset(top_250_genes_exp, select = -c(X1))


library(matrixStats, quietly = TRUE)
nrow(top_250_genes_exp)

#data_mean_filtered <- filterGenes(gse_data_frame, filterTypes = "central", filterCentralType = "median",
#                                 filterCentralPercentile = 0.3)

#nrow(data_mean_filtered)

cor_res <- getCors(inputMat = top_250_genes_exp, design = design_mat)
str(cor_res)

dcPairs_res <- pairwiseDCor(cor_res, compare = c("Normal", "Gastric"))
str(dcPairs_res)

dd_pairs <- dcTopPairs(dcPairs_res, nPairs = 25000, classify = TRUE)

#remove one outlier sample before visualization

#Gets all top dc pairs with CXXC5
cxxc5_pairs <- dd_pairs[dd_pairs$Gene1 %in% cxxc5_probes,]
genes <- c(cxxc5_pairs$Gene2, cxxc5_probes)
#gene_names <- subset(gse54129_ranked_genes[gse54129_ranked_genes$ID %in% genes,],
                     #select=c(ID, Gene.symbol, Gene.title))

#ddcor_res <- ddcorAll(inputMat = top_250_genes_exp, design = design_mat,
                      #compare = c("Normal", "Gastric"), adjust = "none",
                      #heatmapPlot = FALSE, nPerm = 0, splitSet = "222996_s_at")
#head(ddcor_res)

corr <- c("+/+", "+/-", "-/-", "-/+", "0/+", "0/-")
cxxc5_pairs_cleaned <- cxxc5_pairs[cxxc5_pairs$Classes %in% corr,]

#write.csv(genes, file="Probe_to_Symbols.csv")


#a <- cxxc5_pairs_cleaned[, rownames(cxxc5_pairs_cleaned$Gene1) %in% gene_names_from_david]
#gastrics <- GSE66229_expression_matrix[, !names(GSE66229_expression_matrix) %in% 
                                         #normal_samples]

#write.csv(cxxc5_pairs_cleaned, file = "CXXC5_Pairs_Genes_Non-Labeled.csv")



gene1 <- c(CXXC5_Pairs_Genes_Labeled$Gene1)
gene2 <- c(CXXC5_Pairs_Genes_Labeled$Gene2)
corA <- c(CXXC5_Pairs_Genes_Labeled$corA)
corB <- c(CXXC5_Pairs_Genes_Labeled$corB)
classes <- c(CXXC5_Pairs_Genes_Labeled$Classes)

cxxc5_to_table <- data.frame(geneA = gene1, geneB = gene2, cor1 = corA,
                             cor2 = corB, class = classes)

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "229177_at", xlab = "CXXC5", ylab = "C16orf89")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "224516_s_at",
         geneB = "225799_at", xlab = "CXXC5", ylab = "CYTOR")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "224516_s_at",
         geneB = "227198_at", xlab = "CXXC5", ylab = "AFF3")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "224516_s_at",
         geneB = "224654_at", xlab = "CXXC5 ", ylab = "DDX21")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "231015_at", xlab = "CXXC5", ylab = "KLF15")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "224516_s_at",
         geneB = "231015_at", xlab = "CXXC5", ylab = "PLAC9")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "226303_at", xlab = "CXXC5", ylab = "PGM5")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "238075_at", xlab = "CXXC5", ylab = "CHEK1")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "226319_s_at", xlab = "CXXC5", ylab = "ALYREF")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "228802_at", xlab = "CXXC5", ylab = "RBPMS2")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "225915_at", xlab = "CXXC5", ylab = "CAB39L")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "48031_r_at", xlab = "CXXC5", ylab = "FAXDC2")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "224516_s_at",
         geneB = "229019_at", xlab = "CXXC5 ", ylab = "ZNF3858")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "224516_s_at",
         geneB = "229288_at", xlab = "CXXC5 ", ylab = "EPHA7")

plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "224516_s_at",
         geneB = "226192_at", xlab = "CXXC5 ", ylab = "AR")































