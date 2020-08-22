library(DGCA)
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2, quietly= TRUE)
library(formattable)
library(impute, quietly = TRUE)

#Gets rid of quotations in everything
gse33651_ranked_genes <- as.data.frame(sapply(gse33651_ranked_genes, function(x) gsub("\"", "", x)))
#write.table(gse33651_ranked_genes, "gse33651ranked.txt", row.names=FALSE, sep="\t", quote=FALSE)



top_250_ranked_genes <- top_n(gse33651_ranked, -250, P.Value)




n_normal_samples <- 12
n_gastric_samples <- 40

cell_type <- c(rep("Normal", n_normal_samples), rep("Gastric", n_gastric_samples))
design_mat <- model.matrix(~ cell_type + 0)
colnames(design_mat) <- c("Normal", "Gastric")
x <- str(design_mat)

design_mat <- makeDesign(cell_type)
design_mat


top_250_genes_id <- c(top_250_ranked_genes$"ID", "380062")


top_250_genes_exp <- gse_data_frame[row.names(gse_data_frame) %in% top_250_genes_id,]

library(matrixStats, quietly = TRUE)
nrow(top_250_genes_exp)

#data_mean_filtered <- filterGenes(gse_data_frame, filterTypes = "central", filterCentralType = "median",
#                                 filterCentralPercentile = 0.3)

#nrow(data_mean_filtered)

cor_res <- getCors(inputMat = top_250_genes_exp, design = design_mat)
str(cor_res)

dcPairs_res <- pairwiseDCor(cor_res, compare = c("Normal", "Gastric"))
str(dcPairs_res)

dd_pairs <- dcTopPairs(dcPairs_res, nPairs = 250, classify = TRUE)

cxxc5_pairs <- dd_pairs[dd_pairs$Gene1 %in% "380062",]
genes <- c(cxxc5_pairs$Gene2, "380062")
gene_names <- subset(gse33651_ranked[gse33651_ranked$ID %in% genes,],
                     select=c(ID, Gene.symbol, Gene.title))





