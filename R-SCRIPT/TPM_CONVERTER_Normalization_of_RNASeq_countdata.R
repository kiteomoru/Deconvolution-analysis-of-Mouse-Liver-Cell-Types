#Set path
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))



library(biomaRt)
library(tidyverse)
counts= read.delim('/Users/kiteomoru/Downloads/JangaJOB/HTSEQ_COUNTS/all/ALLRAW_COUNTS_JANGA_FURUTA copy.txt', header = T)
counts=as.data.frame(counts)
ncol(counts)
head(counts)
view(counts)
GENES= counts$GeneSymbol
mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene_coords=getBM(attributes=c("mgi_symbol","ensembl_gene_id", "start_position","end_position"), filters="ensembl_gene_id", values=GENES, mart=mouse)
gene_coords$size=gene_coords$end_position - gene_coords$start_position

head(gene_coords)
genes_df= cbind(gene_coords$ensembl_gene_id, gene_coords$size)
genes_df=as.data.frame(genes_df)
head(genes_df)
colnames(genes_df)= c('gene', 'length')

counts2= counts[,2:11]
rownames(counts2)= counts[,1]
#tpm
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

str(counts2)
str(genes_df)
tpms <- apply(counts2, 2, function(x) tpm(x, as.numeric(genes_df$length)))
tpms=as.data.frame(tpms)
view(tpms)
head(tpms)
tpms= as.data.frame(tpms)
write.table(tpms, file = 'TPMSALLFURUTAJANGADATADATA.txt')
