library(Biobase)
library(dplyr)
library(GEOquery)
gse <- getGEO(filename="GSE66229_series_matrix.txt.gz",GSEMatrix = TRUE,getGPL = TRUE) #Retrieve matrix data and store it in R object
show(object = gse) ## To summarize the gse object
x <- exprs(object = gse) #Get expression set from gse object

#write.csv(x = x, file = "GSE66229.expression.matrix.csv", quote = T, row.names = T) #export expression matrix in file (.csv format).
normal_samples <- c(normal_gsm$...1)
norms <- GSE66229_expression_matrix[,c("X1", normal_samples)]
gastrics <- GSE66229_expression_matrix[, !names(GSE66229_expression_matrix) %in% normal_samples]
gse_normal_then_gastric <- merge(norms, gastrics, by="X1")


