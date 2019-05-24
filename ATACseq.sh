#!/bin/bash

# loading additional script to call functions from. Copy for note-taking: /seq/vgb-HAL/swofford/RNAseq/bashFunctionsJasonEdit.sh
. /seq/vgb/jason/scripts/bashFunctions.sh

source /broad/software/scripts/useuse

description="Describes the contents of $(pwd)"
argumentDef=()

# base script used in each project directory, modified to do specific tasks, but also serves as descriptor of project for anyone looking. Following options not directly related to ATACseq.
# The "0" and "no" statements stand for the number of arguments allowed with 
# that option, and the default argument(s) provided

argumentDef+=("-c" printContents "Prints the directory contents" 0 "no")
argumentDef+=("-d" printDescription "Prints a description of the directory's purpose" 0 "no")
argumentDef+=("-fastqc" runFastQC "Runs fastq on all merged fastq files sequencing quality control stats" 0 "no")
argumentDef+=("-compile" compileIndex "Compiles a bowtie2 index for canFam3.1" 0 "no")
argumentDef+=("-revert" revertBams "Reverts the bam files into fastas for alignment.." 0 "no")
argumentDef+=("-align" runAlign "Runs Bowtie2 on the queuing system to align the reads." 0 "no")
argumentDef+=("-i_align" serverAlign "For internal use; code that runs Bowtie2 from SGE" 1 "null")
argumentDef+=("-mergefq" mergeFastQ "Merge multiple lanes of fastqs into one file for alignment" 0 "no")
argumentDef+=("-i_revert" serverRevert "For internal use; code that runs picard on UGER to revert bams to FASTQ" 1 "null")
argumentDef+=("-i_removeM" i_removeM "For internal use; runs samtools/awk on UGER to remove chrM" 1 "null")
argumentDef+=("-removechrM" removeM "remove mitochondrial reads, also reports fraction of reads aligning to chrM" 0 "no")
argumentDef+=("-removePCRdup" dedupPCR "remove PCR duplicats" 0 "no")
argumentDef+=("-i_dedupPCR" i_dedupPCR "For internal use; runs picardtools on UGER to remove PCR duplicates" 1 "null")
argumentDef+=("-removeNonUnique" nonUnique "remove ambiguous (non-unique) alignments" 0 "no")
argumentDef+=("-i_nonUnique" i_nonUnique "For internal use; runs samtools on UGER to remove non-unique alignments" 1 "null")
argumentDef+=("-samToBed" samToBed "convert the alignments to BED intervals, keep PE info, but extending singletons to keep as much info as possible" 0 "no")
argumentDef+=("-i_samToBed" i_samToBed "For internal use; converts bam to bed, extending unpaired reads by average derived from rest of file" 1 "null")
argumentDef+=("-callPeaks" callPeaks "call peaks using MACS2" 0 "no")
argumentDef+=("-cleanup" cleanup "remove fastqs and early step bams from project directory" 0 "no")
argumentDef+=("-runataqv" atacqv "QC ATACseq final results" 0 "no")
argumentDef+=("-runmkarv" mkarv "Collects ataqv results and creates visualiztion html" 0 "no")


# returns the number of arguments passed to the script, not counting the script itself as an argument (stored as $0)
if [ "$#" == "0" ]; then
    argumentsArr=("-h")
else
    # Presumably stores all of the arguments passed to the script in an array named 'argumentsARR' - after checking, this is basically accurate. $@ stores all arguments passed to script as separate, quoted strings. A similar variable "$*" stores all of the arguments as a single string. 
    argumentsArr=("$@")
fi

# handleArguments is a function from /seq/vgb/jason/scripts/bashFunctions.sh
# takes as arguments: 
# description (a description of this script, as defined above.
# argumentDef, which is an array of possible arguments, and their descriptions
# as exemplified above
# argumentsArr, which is an array of input arguments passed to the script,
# as parsed above.  
handleArguments "$description" argumentDef[@] argumentsArr[@]

if [ $printContents == "yes" ]; then
    echo
    echo "dir_contents.sh          This script."
    echo
fi

if [ $printDescription == "yes" ]; then
    echo
    echo "This folder contains files related to assembling and aligning dog ATACseq data"
    echo "Following the best practices from Harvard FAS site, adapting for Broad"
    echo "ref: https://informatics.fas.harvard.edu/atac-seq-guidelines.html"
    echo "Should be broadly applicable to ChIPseq as well, should anyone do that"
    echo "Actual script adapted from RNAseq calling script, may be leftover bits unrelated to ATAC"
    echo "For consistancy and control, will revert bams to fasta and realign, etc."
    echo "General pattern: revert bams, QC sequencing (FASTQC,) align with bowtie2, remove duplicates (Picard,) "
    echo "remove chrM, call peaks with MACS2, QC ATAC (ATACQV,) annotate peaks (chipseeker,)"
    echo "look at motifs (HOMER)"
    echo ""
    echo "Project folder moved from /seq/vgb-HAL/ due to decomissioning. May be some incorrect "
    echo "references left in script from the move."
    echo
fi

picardVersion="/seq/software/picard/1.1178/bin/picard.jar"
gatkVersion="/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-0-gcfedb67/GenomeAnalysisTK.jar"
dogRef="/seq/references/Canis_lupus_familiaris_assembly3/v1/Canis_lupus_familiaris_assembly3.fasta"
ataqvDir="/seq/vgb/swofford/software/ataqv-1.0.0/bin/"


if [ $compileIndex == "yes" ]; then
    reuse UGER
    reuse Bowtie2
    mkdir genome_index
    qsub -cwd -l h_rt=24:00:00 -l h_vmem=20G -b y -V bowtie2-build $dogRef ./genome_index/canFam3.index
fi

if [ $mergeFastQ == "yes" ]; then
    mkdir ./mergedFastqs
    for id in `find ./reverted_fastqs/*.0.fq.gz | sed -e 's:^./reverted_fastqs/::' | sed -e 's:_.*::' | uniq`      
    do
	if [ ! -f ./mergedFastqs/"$id".fwd.fq.gz ];then
	    cat ./reverted_fastqs/"$id"*.0.fq.gz > ./mergedFastqs/"$id".fwd.fq.gz
	    cat ./reverted_fastqs/"$id"*.1.fq.gz > ./mergedFastqs/"$id".rv.fq.gz
	fi
    done
fi

if [ $runFastQC == "yes" ]; then
    cd ./mergedFastqs
    reuse .fastqc-0.11.4
    reuse Picard-Tools
    mkdir ./fastqc
    mkdir ./jobOut
    for file in ./*.fq.gz
    do
	qsub -b y -cwd -V -o ./jobOut/ fastqc $file 
    done
fi

if [ $revertBams == "yes" ]; then
    reuse UGER
    # check to see if subdirectory 'original_bams' exists, if not, create.
    if [ ! -d original_bams ]; then
	mkdir original_bams/
    fi
    cd original_bams/

    ######
    # list of Files upon which to work:
    ######

    for file in /seq/vgb/rawData/ATAC_Seq/productionFirst30_032719/get.broadinstitute.org/pkgs/SN0169398/*.bam

    do
	
	realFile=`readlink -f $file`
	if [ ! -f `basename $realFile` ]; then
	    ln -s $realFile
	fi
    done

    cd ../

    if [ ! -d reverted_fastqs ]; then
	mkdir reverted_fastqs
    fi
    
    if [ ! -d jobOut ]; then
	mkdir jobOut
    fi

    cd reverted_fastqs

    for file in ../original_bams/*
    do
	# checking to see if FASTQs already exist (as .fq.gz) before reverting
	id=`basename $file | sed 's:.bam::'`
	if [ ! -f $id.1.fq.gz ]; then
	    qsub -cwd -o ./jobOut/ -l h_vmem=9G -l h_rt=4:00:00 -j y ../$0 -i_revert $id
	fi
    done
fi

alignmentDir="alignments/"
numCPUs="4"

if [ $runAlign == "yes" ]; then
    reuse Bowtie2
    reuse UGER
    reuse Samtools
    #reuse GridEngine8
    if [ ! -d $alignmentDir ]; then
	mkdir $alignmentDir
    fi
    for file in mergedFastqs/*.fwd.fq.gz
    do 
	id=`basename $file | sed 's:.fwd.fq.gz::'`
	if [ ! -d $alignmentDir/$id ]; then
	    qsub -S /bin/bash -pe smp $numCPUs -binding linear:$numCPUs -V -cwd -o ./jobOut -l h_vmem=10G -l h_rt=24:00:00 -j y $0 -i_align $id 
	fi
    done
fi


if [ $serverRevert != "null" ]; then
    reuse Java-1.8
    id=`echo $serverRevert | sed 's:^[0-9]_::'`
    java -Xmx5G -jar $picardVersion SamToFastq I=../original_bams/$serverRevert.bam FASTQ="$id".0.fq SECOND_END_FASTQ="$id".1.fq UNPAIRED_FASTQ="$id".2.fq >& "$serverRevert".revert.out
    gzip $id.*fq
fi


if [ $serverAlign != "null" ]; then
    reuse Samtools
    mkdir $alignmentDir/$serverAlign/
    cd $alignmentDir/$serverAlign/
    bowtie2 --very-sensitive -p $numCPUs -x /seq/vgb/swofford/ATAC/scripts/genome_index/canFam3.index -X 1000 -1 ../../mergedFastqs/$serverAlign.fwd.fq.gz -2 ../../mergedFastqs/$serverAlign.rv.fq.gz | samtools view -u - | samtools sort - > $serverAlign.bam
    samtools index -b $serverAlign.bam $serverAlign.bam.bai
fi

if [ $removeM == "yes" ]; then
    for file in ./mergedFastqs/*.fwd.fq.gz
    do
	id=$(basename "$file" .fwd.fq.gz)
	if [ ! -f ./$alignmentDir"$id"/"$id"_chrMdel.bam ]; then
	    qsub -V -cwd -l h_vmem=10G -l h_rt=10:00:00 -o ./jobOut/ $0 -i_removeM $id
	fi
    done
fi

if [ $i_removeM != "null" ]; then
    reuse Samtools
    samtools view -h ./$alignmentDir"$i_removeM"/"$i_removeM".bam | awk '{if($3 != "chrM"){print $0}}' | samtools view -Sb - > ./$alignmentDir"$i_removeM"/"$i_removeM"_chrMdel.bam
    samtools index ./$alignmentDir"$i_removeM"/"$i_removeM"_chrMdel.bam
fi

if [ $dedupPCR == "yes" ]; then
    reuse UGER
    for file in ./mergedFastqs/*.fwd.fq.gz
    do
	id=$(basename "$file" .fwd.fq.gz)
	if [ ! -f ./$alignmentDir/"$id"/"$id"_chrMdePCRdup.bam ]; then
	    qsub -V -cwd -o ./jobOut/ -l h_vmem=10G -l h_rt=10:00:00 $0 -i_dedupPCR $id
	fi
    done
fi

if [ $i_dedupPCR != "null" ]; then
    reuse Java-1.8
    java -Xmx8G -jar $picardVersion MarkDuplicates I=./$alignmentDir/"$i_dedupPCR"/"$i_dedupPCR"_chrMdel.bam O=./$alignmentDir/"$i_dedupPCR"/"$i_dedupPCR"_chrMdePCRdup.bam M=./$alignmentDir/"$i_dedupPCR"/dups.txt REMOVE_DUPLICATES=true
fi

if [ $nonUnique == "yes" ]; then
    reuse UGER
    for file in ./mergedFastqs/*.fwd.fq.gz
    do
        id=$(basename "$file" .fwd.fq.gz)
	if [ ! -f ./"$alignmentDir"/"$id"/"$id"_chrMdePCRdupUniqAl.bam ]; then
	    qsub -V -cwd -j y -o ./jobOut/ -l h_vmem=10G -l h_rt=10:00:00 $0 -i_nonUnique $id
	fi
    done
fi

if [ $i_nonUnique != "null" ]; then
    reuse Samtools
    samtools view -b -q 10 ./"$alignmentDir"/"$i_nonUnique"/"$i_nonUnique"_chrMdePCRdup.bam > ./"$alignmentDir"/"$i_nonUnique"/"$i_nonUnique"_chrMdePCRdupUniqAl.bam
    samtools index ./"$alignmentDir"/"$i_nonUnique"/"$i_nonUnique"_chrMdePCRdupUniqAl.bam
fi

if [ $samToBed == "yes" ]; then
    for file in ./mergedFastqs/*.fwd.fq.gz
    do
        id=$(basename "$file" .fwd.fq.gz)
	if [ ! -f ./"$alignmentDir"/"$id"/"$id"_chrMdePCRdupUniqAl.bed ]; then
	    qsub -V -cwd -j y -o ./jobOut/ -l h_vmem=10G -l h_rt=10:00:00 $0 -i_samToBed $id
	fi
    done
fi

if [ $i_samToBed != "null" ]; then
    reuse Samtools
    reuse Anaconda3
    samtools view -h ./"$alignmentDir"/"$i_samToBed"/"$i_samToBed"_chrMdePCRdupUniqAl.bam | python SAMtoBED.py -i - -o ./"$alignmentDir"/"$i_samToBed"/"$i_samToBed"_chrMdePCRdupUniqAl.bed -x -v
fi 

if [ $callPeaks == "yes" ]; then
    reuse UGER
    reuse Anaconda
    source activate /seq/vgb/swofford/conda/my2env/
    for file in ./mergedFastqs/*.fwd.fq.gz
    do
	id=$(basename "$file" .fwd.fq.gz)
	if [ ! -d ./"$alignmentDir"/"$id"/peaks ]; then
	    qsub -V -cwd -l h_rt=5:00:00 -l h_vmem=5G -j y -o ./jobOut/ -b y macs2 callpeak -t ./"$alignmentDir"/"$id"/"$id"_chrMdePCRdupUniqAl.bed -f BEDPE -n $id -g 2.5e9 --keep-dup all --outdir ./"$alignmentDir"/"$id"/peaks/
	fi
    done
fi

if [ $atacqv == "yes" ]; then
    reuse UGER
    for file in ./mergedFastqs/*.fwd.fq.gz
    do
        id=$(basename "$file" .fwd.fq.gz)
	if [ ! -d ./"$alignmentDir"/"$id"/peaks/QC ]; then
	    mkdir ./"$alignmentDir"/"$id"/peaks/QC
	    qsub -cwd -V -l h_rt=10:00:00 -l h_vmem=5G -j y -o ./jobOut/ -b y "$ataqvDir"ataqv --autosomal-reference-file /seq/vgb/swofford/ref/chr.txt --mitochondrial-reference-name chrM --tss-file /seq/vgb/swofford/ref/TSSFullUniqueGene.bed --metrics-file ./"$alignmentDir"/"$id"/peaks/QC/"$id".ataqv.json --ignore-read-groups --name $id  --peak-file ./"$alignmentDir"/"$id"/peaks/"$id"_peaks.narrowPeak Dog ./"$alignmentDir"/"$id"/"$id".bam
	fi
    done
fi

if [ $mkarv == "yes" ]; then
    reuse UGER
    reuse Anaconda
    source activate /seq/vgb/swofford/conda/my2env/

    if [ ! -d ./QCVisualization ]; then
	mkdir ./QCVisualiztion
    fi
    
    qsub -V -cwd -l h_vmem=50G -l h_rt=10:00:00 -b y -j y -o ./jobOut/ "$ataqvDir"mkarv --force ./QCVisualization ./"$alignmentDir"/*/*/*/*.json
fi
	



