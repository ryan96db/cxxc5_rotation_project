library(DGCA, quietly = TRUE)
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2, quietly= TRUE)
library(formattable)
library(data.table)


ranked_genes <- fread("gse19826_ranked_genes.txt")
top_250 <- head(ranked_genes, 250)



n_normal_samples <- 15
n_gastric_samples <- 12

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

top_250_genes_exp <- as.data.frame(normal_then_gastric[normal_then_gastric$X1 %in% top_250_genes_id,])
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


cxxc5_pairs <- dd_pairs[dd_pairs$Gene1 %in% cxxc5_probes,]
genes <- c(cxxc5_pairs$Gene2, cxxc5_probes)

corr <- c("+/+", "+/-", "-/-", "-/+", "0/+", "0/-")
cxxc5_pairs_cleaned <- cxxc5_pairs[cxxc5_pairs$Classes %in% corr,]

#write.csv(genes, file="Probe_to_Symbols.csv")
gene_names_from_david <- data.frame(david)

a <- merge(cxxc5_pairs_cleaned, gene_names_from_david, by="Gene2")

gene1 <- c(a$Gene1)
gene2 <- c(a$Gene2)
corA <- c(a$corA)
corB <- c(a$corB)
classes <- c(a$Classes)
corB_p <- c(a$corB_pVal)

cxxc5_to_table <- data.frame(geneA = gene1, geneB = gene2, cor1 = corA,
                             cor2 = corB, class=classes, corB_p = corB_pVal)



