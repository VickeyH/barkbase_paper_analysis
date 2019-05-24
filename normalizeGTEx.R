library(edgeR)

#Read in the raw counts (be sure set row.names if list of genes)
counts <- read.delim("/seq/vgb/barkbase/human_comparison/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct", row.names="Name",sep="\t")
counts <- subset(counts, select=-c(Description))

#Make a DGEList object
dge <- DGEList(counts=counts)

#Filter genes not expressed at cpm >1 in >=2 samples
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors
dge <- calcNormFactors(dge, method="TMM")

#Calculate cpms
cpm <- cpm(dge)

#Output table (if needed)
write.table(cpm, file="GTEx_TMM_norm_gene_counts.txt", sep="\t", col.names=NA, quote=FALSE)
