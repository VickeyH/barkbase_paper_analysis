# Read in human count data
counts <- read.table("GTEx_TMM_norm_gene_counts.txt", sep="\t", header=T, stringsAsFactors=F, row.names=1)
attributes <- read.table("GTEx_v7_Annotations_SampleAttributesDS.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
# For each tissue that we are interested in
manifest <- read.table("/seq/vgb/barkbase/embryo_DGE/sample.manifest.cleaned", header=F, stringsAsFactors=F, sep="\t")
manifest <- manifest[,c(1,2,3,10)]
colnames(manifest) <- c("publicname", "sample", "sampleid", "tissue")
tissues <- unique(manifest$tissue)
tissues.human <- unique(attributes$SMTSD)
tissues.map <- read.table("tissuemapping.txt", sep="|", header=T, stringsAsFactors=F)
counts.medians <- data.frame(ens_id = rownames(counts))
for(tissue in tissues) {
  counts.medians <- cbind(counts.medians, rep(NA, nrow(counts.medians)))
}
colnames(counts.medians) <- c("gene_id", tissues)
for (tissue in tissues) {
  # Get the list of human samples from that tissue
  tissue.human <- tissues.map$human[tissues.map$dog == tissue]
  if (is.na(tissue.human)) { next }
  samples.human <- attributes$SAMPID[attributes$SMTSD == tissue.human]
  samples.human <- gsub("-", "\\.", samples.human, perl=TRUE)
  # For each gene (V1) in counts, get median of counts for each sample from this tissue
  counts.tissue <- sapply(1:nrow(counts), function(row) {
    median(as.numeric(counts[row, colnames(counts) %in% samples.human]))
  } )
  counts.medians[,tissue] <- counts.tissue
}
write.table(counts.medians, file="human-count-medians-output.txt", row.names=F, colnames=T, quote=F)
