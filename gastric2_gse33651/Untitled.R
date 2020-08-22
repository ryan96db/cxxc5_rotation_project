library(Biobase)
library(dplyr)
library(GEOquery)
gse54129 <- getGEO(filename = 'GSE33651_family.soft.gz')
#x <- Table(GSMList(gse54129)[[1]])[1:5,]
gsm_names <- names(GSMList(gse54129))
len <- length(gsm_names)

gse_id_ref <- c(Table(GSMList(gse54129)[[1]])$ID_REF)

gse_data_frame <- data.frame(matrix(ncol = len, nrow = length(gse_id_ref)))

gse_names <- c(gsm_names)

colnames(gse_data_frame) <- gse_names
rownames(gse_data_frame) <- gse_id_ref

