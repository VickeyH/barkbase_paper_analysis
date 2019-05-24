library(GenomicFeatures)
library(ChIPseeker)

#build reference database for annotation
TxDb.canFam3.EnsemblGene <- makeTxDbFromUCSC(genome = 'canFam3', tablename = "ensGene")

#load list of peak files (from MACS2)
tissuepeaks <- dir("/seq/vgb/swofford/ATAC/peaks/Tissue/", pattern = "*.narrowPeak", full.names = TRUE)

#define promotor windows
promoter <- getPromoters(TxDb=TxDb.canFam3.EnsemblGene, upstream=3000, downstream=3000)

#create tag matrix for promoter windows
tagMatrixList <- lapply(tissuepeaks, getTagMatrix, windows=promoter)
names(tagMatrixList) <- gsub('.narrowPeak','',basename(tissuepeaks))

#plot count frequency over promotor window range
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), names=names(tagMatrixList))

#annotate peaks
tissueAnno <- lapply(tissuepeaks, annotatePeak, TxDb=TxDb.canFam3.EnsemblGene, tssRegion=c(-3000,3000), verbose=FALSE)
names(tissueAnno) <- gsub('.narrowPeak','',basename(tissuepeaks))

#plot Annotation type
plotAnnoBar(tissueAnno)

#plot distance to Transcription Start Site
plotDistToTSS(tissueAnno)

#load specific peak file for annotation plotting
PancreasPeak <- "/seq/vgb/swofford/ATAC/peaks/Tissue/pancreas.narrowPeak"
PancreasAnno <- annotatePeak(PancreasPeak, tssRegion=c(-3000,3000), TxDb=TxDb.canFam3.EnsemblGene, verbose=FALSE)

#upset plot
upsetPlot(PancreasAnno)