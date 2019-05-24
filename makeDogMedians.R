# Read in dog count data
counts <- read.table("TMM_norm_gene_counts.txt", sep="\t", header=T, stringsAsFactors=F, row.names=1)
attributes <- read.table("/seq/vgb/barkbase/embryo_DGE/sample.manifest.cleaned", sep="\t", header=F, stringsAsFactors=F)
attributes <- attributes[,c(1,2,3,10)]
colnames(Attributes) <- c("publicname", "sample", "sampleid", "tissue")

# For each tissue that we are interested in
#manifest <- read.table("/seq/vgb/barkbase/embryo_DGE/sample.manifest.cleaned", header=F, stringsAsFactors=F, sep="\t")
#manifest <- manifest[,c(1,2,3,10)]
#colnames(manifest) <- c("publicname", "sample", "sampleid", "tissue")

#tissues <- unique(manifest$tissue)
#tissues.human <- unique(attributes$SMTSD)
#tissues.map <- read.table("tissuemapping.txt", sep="|", header=T, stringsAsFactors=F)

tissues <- unique(attributes$tissue)

counts.medians <- data.frame(ens_id = rownames(counts))
for(tissue in tissues) {
  counts.medians <- cbind(counts.medians, rep(NA, nrow(counts.medians)))
}
colnames(counts.medians) <- c("gene_id", tissues)

for (tissue in tissues) {
  # Get the list of human samples from that tissue
  #tissue.human <- tissues.map$human[tissues.map$dog == tissue]
  #if (is.na(tissue.human)) { next }

  samples.dog <- attributes$sample[attributes$tissue == tissue]
  #samples.human <- attributes$SAMPID[attributes$SMTSD == tissue.human]
  #samples.human <- gsub("-", "\\.", samples.human, perl=TRUE)

  # For each gene (V1) in counts, get median of counts for each sample from this tissue
  counts.tissue <- sapply(1:nrow(counts), function(row) {
    median(as.numeric(counts[row, colnames(counts) %in% samples.dog]))
  } )
  counts.medians[,tissue] <- counts.tissue
}

write.table(counts.medians, file="dog-count-medians-output.txt", row.names=F, col.names=T, quote=F)
