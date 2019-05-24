library(ggplot2)

RNACounts <- read.table('/Volumes/seq_vgb/barkbase/counts/TMM/TMM_norm_transcript_counts.txt')
#make a true/false list based on tissue names
isPancreas <- grepl('ancreas', names(RNACounts))
isSalivary <- grepl('aliva', names(RNACounts))

#make new dataframe with only that tissue
pancreasCounts <- RNACounts[,isPancreas]
salivaryCounts <- RNACounts[,isSalivary]

novelTranscripts <- read.table('/Volumes/seq_vgb/barkbase/exome_comparison/list_unmached_in_ensembl_hoeppner_multi_exon.txt')
cpmCutoff <- 1
#dumb way of getting the novel Transcripts into a list
novelTranscripts <- novelTranscripts$V1

#make new data frame that only has novel genes
pancreasNovel <- subset(pancreasCounts, rownames(pancreasCounts) %in% novelTranscripts)
salivaryNovel <- subset(salivaryCounts, rownames(salivaryCounts) %in% novelTranscripts)

#get the mean cpm (change the 5 here if tissue has a different number of dogs)
pancreasNovel$mean <- rowMeans(pancreasNovel[,1:5])
salivaryNovel$mean <- rowMeans(salivaryNovel[,1:5])

#count dogs with expression over 0, again, change the 5 here if needed
pancreasNovel$dogCount <- apply(pancreasNovel[,1:5] > 0, 1, sum)
salivaryNovel$dogCount <- apply(salivaryNovel[,1:5] > 0, 1, sum)

#new column with pass/fail for mean cpm >1 and more than one dog
pancreasNovel$pass <- ifelse(pancreasNovel$mean > cpmCutoff & pancreasNovel$dogCount > 1,'pass','fail')
salivaryNovel$pass <- ifelse(salivaryNovel$mean > cpmCutoff & salivaryNovel$dogCount > 1, 'pass','fail')

#new df with only passed genes
pancreasNovelKeep <- pancreasNovel[pancreasNovel$pass == 'pass',]
salivaryNovelKeep <- salivaryNovel[salivaryNovel$pass == 'pass',]

#list of novel transcripts from row names of passed genes
pancreasNovelTranscripts <- row.names(pancreasNovelKeep)
salivaryNovelTranscripts <- row.names(salivaryNovelKeep)

#get just the mstgr.blah13 bit of the gene name
pancreasNovelGenes <- unique(sub('([A-Z]+.[0-9]+).[0-9]+','\\1',pancreasNovelTranscripts))
salivaryNovelGenes <- unique(sub('([A-Z]+.[0-9]+).[0-9]+','\\1',salivaryNovelTranscripts))

#novel genes with zero tissue expression
nonPancreasNovel <- pancreasNovel[pancreasNovel$mean == 0,]
nonPancreasNovelTranscripts <- rownames(nonPancreasNovel)

nonPancreasNovelGenes <- unique(sub('([A-Z]+.[0-9]+).[0-9]+','\\1',nonPancreasNovelTranscripts))

stringtieGTF <- read.table('/Volumes/seq_vgb-hal/swofford/RNAseq/assemblies/stringtie/stringMerged/stringtie_Merged.gtf', sep = '\t')

#dump everything but the Transcript lines
stringtieTranscripts <- stringtieGTF[stringtieGTF$V3 == 'transcript',]
stringtieTranscripts$gene <- sub('^gene_id\\s(MSTRG.[0-9]+);\\s*.*','\\1',stringtieTranscripts$V9)

#keep only the novel genes from tissue
stringtieKeep <- subset(stringtieTranscripts, stringtieTranscripts$gene %in% pancreasNovelGenes)
stringtieKeepSalivary <- subset(stringtieTranscripts, stringtieTranscripts$gene %in% salivaryNovelGenes)

#get columns for bed file
stringBed <- stringtieKeep[c('V1','V4','V5','gene','V8','V7')]
stringSalivaryBed <- stringtieKeepSalivary[c('V1','V4','V5','gene','V8','V7')]
write.table(stringBed, file = '/Volumes/seq_vgb/swofford/R24/novelGenesPancreas.bed', col.names = F, row.names = F, quote = F, sep = '\t')
write.table(stringSalivaryBed, file = '/Volumes/seq_vgb/swofford/R24/novelGenesSalivary.bed', col.names = F, row.names = F, quote = F, sep = '\t')


#bash and bedtools to sort, merge, get distance
#sort -k1,1 -k2,2n bedFile > sortedBedFile
#bedtools merge -i sortedBedFile > mergedBedFile
#bedtools closest -a mergedBedFile -b mergedSortedPeakFile -D ref > closest.out
#bedtools merge -c 22 -o distinct -i closest.out > closestMerged.out 



distNovelPancGenes <- read.table('/Volumes/seq_vgb/swofford/R24/novelTranscripts/novelPancreasDistance.out')
distSalivary <- read.table('/Volumes/seq_vgb/swofford/R24/novelTranscripts/novelSalivaryPancDistanceMerge.out')
distNovelNonPancGenes <- read.table('/Volumes/seq_vgb/swofford/R24/novelTranscripts/novelNonPancreasDistanceMerged.out')

#remove non-numeric
#distNovelNonPancGenes <- distNovelNonPancGenes[distNovelNonPancGenes$V11 !='.',]
distNovelNonPancGenes <- distNovelNonPancGenes[as.numeric(distNovelNonPancGenes$V11) >= 0,]
distNovelPancGenes <- distNovelPancGenes[as.numeric(distNovelPancGenes$V11) >= 0,]
#distNovelPancGenes <- distNovelPancGenes[distNovelPancGenes$V11 !='.',]

#determine whether upstream or down, add directionality to distance
distNovelNonPancGenes$side <- ifelse(distNovelNonPancGenes$V8 < distNovelNonPancGenes$V2,-1,1)
distNovelNonPancGenes$dist <- as.numeric(distNovelNonPancGenes$V11) * as.numeric(distNovelNonPancGenes$side)
distNovelPancGenes$side <- ifelse(distNovelPancGenes$V8 < distNovelPancGenes$V2,-1,1)
distNovelPancGenes$dist <- as.numeric(distNovelPancGenes$V11) * as.numeric(distNovelPancGenes$side)
distSalivary$side <- ifelse(distSalivary$V8 < distSalivary$V2,-1,1)
distSalivary$dist <- distSalivary$V11 * distSalivary$side

#make labeled lists of distances
pancDist <- distNovelPancGenes$dist
pancLabs <- rep('Pancreas',length(pancDist))

nonPancDist <- distNovelNonPancGenes$dist
nonPancLabs <- rep('Non-Pancreas', length(nonPancDist))


salDist <- distSalivary$V11
salLabs <- rep('Salivary', length(salDist))

#combine distances and labels
combName <- c(pancLabs,nonPancLabs)
combDist <- c(pancDist,nonPancDist)

#merge distances and labels into dataframe for plotting
forPlot <- data.frame(Name=as.factor(combName),Distance=combDist)

#plot with ggplot2
myPlot <- ggplot(forPlot, aes(x=Distance, y=Name))  + geom_bin2d(aes(fill = ..density..), bins=100)
myPlot + scale_fill_continuous(low="yellow", high="red") + theme_classic()
