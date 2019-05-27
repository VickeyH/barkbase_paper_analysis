# barkbase_paper_analysis
Analysis scripts relating to the Barkbase paper.

NOTE: These files are presented as-is, for informational purposes.
In other words, while they accurately represent the code that was run to produce the paper analyses
they will require significant rewrites to function outside of our original development environment
as they contain hard coded path information and references to environment-specific packages.

## File descriptions, loosely grouped by analysis:

### RNAseq processing (alignment, assembly, abundance estimation)
- RNAseq alignment and assembly
  - RNAseq.sh
- Calculating raw read counts from StringTie output (using the prepDE.py script from the StringTie online manual)
  - calc_read_count.txt
- Calculating TMM-normalized CPM values in edgeR (and filtering transcripts not expressed at at least 1 CPM in 2 or more dogs)
  - CPM_cutoff_1.txt
- Calculating TMM-normalized CPM values in edgeR (and filtering transcripts not expressed at at least 0.16 CPM in 2 or more dogs)
  - CPM_cutoff_0.16.txt

### ATACseq

- Alignment and peak calling
  - ATACseq.sh
- ATACseqQC plotting (figure 8.)
  - ATACQCplotting.R
- Novel RNA ATACseq proximity (figure 9.)
  - novelTranscriptsToATAC.R

### Human / Dog expression comparison

- Normalization of human read counts
  - normalizeGTEx.R
- Production of human tissue medians
  - makeHumanMedians.R
- Production of dog tissue medians
  - makeDogMedians.R
- Production of spearman matrix
  - produceComparableRankedFiles.py
  
### Ungrouped

- GoSeq enrichment analysis
  - GoSeq_enrichment.txt
- Calculating relatedness using Plink’s pi-hat score
  - relatedness_plink.txt
