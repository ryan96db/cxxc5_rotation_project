library(DGCA)
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2, quietly= TRUE)
library(formattable)

str(Top250Genes, list.len = 5)


n_normal_samples <- 21
n_gastric_samples <- 111

cell_type <- c(rep("Normal", n_normal_samples), rep("Gastric", n_gastric_samples))
design_mat <- model.matrix(~ cell_type + 0)
colnames(design_mat) <- c("Normal", "Gastric")
x <- str(design_mat)

design_mat <- makeDesign(cell_type)
design_mat

top_250_genes_id <- c(Top250Genes$ID, "222996_s_at")

top_250_genes_exp <- gse_data_frame[row.names(gse_data_frame) %in% top_250_genes_id,]

#top_500_genes_exp


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
head(dd_pairs)

#remove one outlier sample before visualization

#Gets all top dc pairs with CXXC5
cxxc5_pairs <- dd_pairs[dd_pairs$Gene1 %in% "222996_s_at",]
genes <- c(cxxc5_pairs$Gene2, "222996_s_at")
gene_names <- subset(gse54129_ranked_genes[gse54129_ranked_genes$ID %in% genes,],
                     select=c(ID, Gene.symbol, Gene.title))

#Adds Gene symbols to CXXC5 gene pairs
cxxc5_pairs_with_gene_symbols <- left_join(cxxc5_pairs, gene_names, by = c("Gene2" = "ID"))
for (inde in 1: length(cxxc5_pairs_with_gene_symbols$Gene1))
{
  cxxc5_pairs_with_gene_symbols$Gene1 <- "CXXC5"
  cxxc5_pairs_with_gene_symbols$Gene2[inde] <- cxxc5_pairs_with_gene_symbols$Gene.symbol[inde]
}

cleaned_up_DGCA_data <- subset(cxxc5_pairs_with_gene_symbols, cxxc5_pairs_with_gene_symbols$Gene2 != "")
cleaned_up_DGCA_data <- cleaned_up_DGCA_data[, -c(10:11)]

missing <- list()

for (indx in 1:length(cxxc5_pairs_with_gene_symbols$Gene2))
{
    if (cxxc5_pairs_with_gene_symbols$Gene.title[indx] == '')
    {
      missing <- append(missing, cxxc5_pairs_with_gene_symbols$Gene2[indx])
      
    }
  
}
missing
#Writes probe ids with no gene symbol to txt file for submission to David
write.table(missing, file = "genes_with_missing_symbol.txt", sep = "\n")

write.csv(cxxc5_pairs_with_gene_symbols, file = "CXXC5_Gene_Pairs.csv")



plotCors(inputMat = top_250_genes_exp, design = design_mat,
         compare = c("Normal", "Gastric"), geneA = "222996_s_at",
         geneB = "225820_at", xlab = "CXXC5", ylab = "JADE1")


