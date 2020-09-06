library(Biobase)
library(dplyr)
library(GEOquery)
gse19826 <- getGEO(filename = 'GSE19826_family.soft.gz')
#x <- Table(GSMList(gse54129)[[1]])[1:5,]
gsm_names <- names(GSMList(gse19826))
len <- length(gsm_names)

gse_id_ref <- c(Table(GSMList(gse19826)[[1]])$ID_REF)

gse_data_frame <- data.frame(matrix(ncol = len, nrow = length(gse_id_ref)))

gse_names <- c(gsm_names)

colnames(gse_data_frame) <- gse_names
rownames(gse_data_frame) <- gse_id_ref


for (ind in 1:len)
{
  sample <- Table(GSMList(gse19826)[[ind]])
  
  sample_name <- gsm_names[ind]
  sample_values <- c(sample$VALUE)
  gse_data_frame[, sample_name] <- sample_values
  
  
}



norm_ <- c(normal$Sample)
norms <- gse_data_frame[which(colnames(gse_data_frame) %in% norm_)]

gastr_ <- c(gastric$Sample)
gastr <- gse_data_frame[which(colnames(gse_data_frame) %in% gastr_)]

write.csv(norms, file="normal_samples_frame.csv")
write.csv(gastr, file="gastric_samples_frame.csv")






