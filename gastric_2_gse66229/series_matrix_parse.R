library(data.table)
library(foreach)

args <- commandArgs(TRUE)

#### Input files ####
# Series Matrix File
SERIES_FILE <- args[1]
# Platform data table obtained from GEO.
PLATFORM_FILE <- args[2]

#### Parameters ####
# Column containing gene IDs in platform file #
# Generally "ENTREZ_GENE_ID" or "GENE"
GENE_COLUMN <- if(!is.na(args[3])) args[3] else "ENTREZ_GENE_ID"

#### Output files ####
EXPRESSION_OUTPUT_FILE <- if(!is.na(args[4])) args[4] else "expression.csv"
ANNOTATION_OUTPUT_FILE <- if(!is.na(args[5])) args[5] else "annotation.csv"

# Read characteristics
con <- file(SERIES_FILE, "r")
characteristics <- c()
while(TRUE) {
  line <- readLines(con, n=1)
  if(length(line) == 0) {
    break
  } else if(startsWith(line, "!Sample_title")) {
    titles <- unlist(strsplit(line, "\t"))[-1]
    titles <- gsub("\\\"", "", titles)
  } else if(startsWith(line, "!Sample_characteristics")) {
    characteristics <- c(characteristics, line)
  } else if(startsWith(line, "!Sample_geo_accession")) {
    accession <- unlist(strsplit(line, "\t"))[-1]
    accession <- gsub("\\\"", "", accession)
  }
}
close(con)

# Parse characteristics
anno <- data.frame(lapply(characteristics, function(x) {
  values <- unlist(strsplit(x, "\t"))[-1]
  values <- gsub("\\\"", "", values)
  parts <- strsplit(values, ": ")
  
  name <- parts[[1]][[1]]
  values <- sapply(parts, function(x) x[2])
  
  out <- list()
  out[[name]] <- values
  return(out)
}))

anno <- data.table(sample=accession, title=titles, anno)

# Read probe-level expression data
D <- fread(SERIES_FILE, header=TRUE, skip="\"ID_REF\"", fill=TRUE, na.strings=c("","NA","null"))
D <- D[1:(nrow(D)-1),] # remove table end marker

# Read platform data table
ref <- read.table(PLATFORM_FILE, header=TRUE, sep="\t", quote="", comment.char="#", fill=TRUE)[,c("ID",GENE_COLUMN)]
colnames(ref) <- c("ID_REF","entrez")
ref$entrez <- as.character(ref$entrez)
ref <- subset(ref, !is.na(entrez) & entrez != "")

entrez_split <- strsplit(ref$entrez, " /// ")
ref <- data.frame(
  ID_REF=rep(ref$ID_REF, sapply(entrez_split, length)),
  entrez=unlist(entrez_split)
)

# Merge tables to map Entrez genes ids
m <- data.table(merge(ref, D, all=FALSE)[,-1])

# Aggregate duplicate genes by median
m <- m[, lapply(.SD, median(na.rm=TRUE)), by=entrez]

# Extract gene and sample names
genes <- m$entrez
samples <- colnames(m)[-1]

# Transpose matrix, add sample column and set gene names as column names
m <- data.table(samples, transpose(m[,-1]))
colnames(m) <- c("sample", as.character(genes))

# Write results to separate expression and annotation files
fwrite(m, file=EXPRESSION_OUTPUT_FILE, sep=",")
fwrite(anno, file=ANNOTATION_OUTPUT_FILE, sep=",")