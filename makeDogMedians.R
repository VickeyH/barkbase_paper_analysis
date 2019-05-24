# Read in dog count data
counts <- read.table("TMM_norm_gene_counts.txt", sep="\t", header=T, stringsAsFactors=F, row.names=1)
attributes <- read.table("/seq/vgb/barkbase/embryo_DGE/sample.manifest.cleaned", sep="\t", header=F, stringsAsFactors=F)
attributes <- attributes[,c(1,2,3,10)]
colnames(Attributes) <- c("publicname", "sample", "sampleid", "tissue")

tissues <- unique(attributes$tissue)

counts.medians <- data.frame(ens_id = rownames(counts))
for(tissue in tissues) {
  counts.medians <- cbind(counts.medians, rep(NA, nrow(counts.medians)))
}
colnames(counts.medians) <- c("gene_id", tissues)

for (tissue in tissues) {
  samples.dog <- attributes$sample[attributes$tissue == tissue]

  # For each gene (V1) in counts, get median of counts for each sample from this tissue
  counts.tissue <- sapply(1:nrow(counts), function(row) {
    median(as.numeric(counts[row, colnames(counts) %in% samples.dog]))
  } )
  counts.medians[,tissue] <- counts.tissue
}

write.table(counts.medians, file="dog-count-medians-output.txt", row.names=F, col.names=T, quote=F)
