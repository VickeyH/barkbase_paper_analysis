#!/bin/bash

# loading additional script to call functions from. Copy for note-taking: /seq/vgb-HAL/swofford/RNAseq/bashFunctionsJasonEdit.sh
. /seq/vgb/jason/scripts/bashFunctions.sh
source /broad/software/scripts/useuse

description="Describes the contents of $(pwd)"
argumentDef=()


argumentDef+=("-c" printContents "Prints the directory contents" 0 "no")
argumentDef+=("-d" printDescription "Prints a description of the directory's purpose" 0 "no")
argumentDef+=("-compile" compileIndex "Compiles a hisat2 index for canFam3.1" 0 "no")
argumentDef+=("-revert" revertBams "Reverts the bam files into fastas for alignment.." 0 "no")
argumentDef+=("-align" runAlign "Runs HISAT2 on the queuing system to align the reads." 0 "no")
argumentDef+=("-i_align" serverAlign "For internal use; code that runs HISAT2 from SGE" 1 "null")
argumentDef+=("-i_revert" serverRevert "For internal use; code that runs picard on UGER to revert bams to FASTQ" 1 "null")
argumentDef+=("-stringtie" runStringtie "Runs stringtie to assemble the separate samples" 0 "no")
argumentDef+=("-stringMerge" runStringMerge "Runs stringmerge to produce one transcript file" 0 "no")
argumentDef+=("-stringQuant" runStringQuant "Estimate abundance using stringtie." 0 "no")
argumentDef+=("-assemblyFASTA" assemblyFASTA "generates FASTA for all assemblies built by stringtie." 0 "no")
argumentDef+=("-fastQC" runFastQC "runs fastqc on all of the hisat aligned bams" 0 "no")
argumentDef+=("-markDuplicates" dedupPCR "runs mark duplicates on samples" 0 "no")
argumentDef+=("-i_dedupPCR" i_dedupPCR "internal, submits mark duplicates to UGER" 1 "null")

if [ "$#" == "0" ]; then
    argumentsArr=("-h")
else
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
    echo "This folder contains files related to assembling and aligning dog RNAseq data"
    echo "General process from bam as starting point is: compile hisat2 ref, revert to fastq, " 
    echo "align with hisat2, stringtie to assemble, stringmerge, and finally stringQuant."
    echo "fasQC can be run to check quality of sequencing, "
    echo "assemblyFASTA will return sequences for all transcripts."
    echo
fi

picardVersion="/seq/software/picard/1.1178/bin/picard.jar"
gatkVersion="/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-0-gcfedb67/GenomeAnalysisTK.jar"
hisatVersion="/seq/vgb/software/hisat2/hisat2-2.1.0/"
dogRef="/seq/references/Canis_lupus_familiaris_assembly3/v1/Canis_lupus_familiaris_assembly3.fasta"
stringtieSoftwareDir="/seq/vgb/software/stringtie/stringtie-1.3.4d.Linux_x86_64/"

if [ $compileIndex == "yes" ]; then
    reuse Python-3.4
    mkdir genome_index
    $hisatVersion/hisat2-build $dogRef genome_index/canFam3.1.plus_MT
fi

if [ $revertBams == "yes" ]; then
    reuse UGER

    if [ ! -d original_bams ]; then
	mkdir original_bams/
    fi
    cd original_bams/

    ######
    # list of Files upon which to work:
    ######


    for file in /seq/picard_aggregation/G144040/*/current/*.bam /seq/picard_aggregation/G130164/*/current/*.bam /seq/picard_aggregation/G145546/*/current/*.bam /seq/picard_aggregation/G144976/*/current/*.bam /seq/picard_aggregation/G110857/*/current/*.bam /seq/picard_aggregation/G124847/*/current/*.bam

    #for file in /seq/picard_aggregation/G110857/*/current/*.bam /seq/picard_aggregation/G124847/*/current/*.bam
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
    cd reverted_fastqs
    for file in ../original_bams/*
    do
	# checking to see if FASTQs already exist (as .fq.gz) before reverting
	id=`basename $file | sed 's:.bam::'`
	echo $id
	if [ ! -f $id.1.fq.gz ]; then
	    qsub -cwd -N $id.revert -l h_vmem=9G -l h_rt=4:00:00 -j y /seq/vgb-HAL/swofford/RNAseq/dir_contentsEDIT.sh -i_revert $id


	fi
    done
fi

alignmentDir="alignments/"
numCPUs=4
if [ $runAlign == "yes" ]; then
    reuse Samtools
    reuse UGER
    if [ ! -d $alignmentDir ]; then
	mkdir $alignmentDir
    fi
    for file in reverted_fastqs/*.0.fq.gz
    do
	id=`basename $file | sed 's:.0.fq.gz::'`
	if [ ! -d $alignmentDir/$id ]; then
	    
	    qsub -S /bin/bash -pe smp $numCPUs -binding linear:$numCPUs -V -cwd -N $id.align -l h_vmem=7G -l os=RedHat7 -j y -o $alignmentDir/$id.out $0 -i_align $id
	    
	fi
    done
fi

if [ $serverRevert != "null" ]; then
    reuse Java-1.8
    java -Xmx5G -jar $picardVersion SamToFastq I=../original_bams/$serverRevert.bam FASTQ="$serverRevert".0.fq SECOND_END_FASTQ="$serverRevert".1.fq UNPAIRED_FASTQ="$serverRevert".2.fq >& "$serverRevert".revert.out
    gzip $serverRevert.*fq
fi

if [ $serverAlign != "null" ]; then
    reuse Samtools
    mkdir $alignmentDir/$serverAlign/
    cd $alignmentDir/$serverAlign/
    $hisatVersion/hisat2 -p $numCPUs --dta --rna-strandness "RF" -x ../../genome_index/canFam3.1.plus_MT -1 ../../reverted_fastqs/$serverAlign.0.fq.gz -2 ../../reverted_fastqs/$serverAlign.1.fq.gz | samtools sort -@ $numCPUs -m 500M -T /seq/vgb/swofford/temp/ -o $serverAlign.sorted.bam
    samtools index -b $serverAlign.sorted.bam $serverAlign.sorted.bam.bai
fi

assemblyDir="assemblies/"
stringtieDir="$assemblyDir/stringtie/"


if [ $runStringtie == "yes" ]; then
    if [ ! -d $assemblyDir ]; then
        mkdir $assemblyDir
    fi

    if [ ! -d $stringtieDir ]; then
        mkdir $stringtieDir
    fi


    reuse UGER
    for file in $alignmentDir/*/*bam
    do
        id=`basename $file | sed 's:.sorted.bam::'`
        if [ ! -d $stringtieDir/$id/ ]; then
            mkdir $stringtieDir/$id/

	     qsub -cwd -b y -V -N $id.stringtie -l h_vmem=10G -j y -o $stringtieDir/$id/$id.stringtie.out $stringtieSoftwareDir/stringtie -o $stringtieDir/$id/$id"_stringtie" -x chrM -G /seq/vgb/swofford/ref/Canis_familiaris.CanFam3.1.95.chrRenamed.gtf -B $file

        fi
    done
fi 



function tarOldVersion()
{
    targetFolder=$1
    archiveFolder=$2
    if [ -d $targetFolder ]; then
        if [ ! -d $archiveFolder ]; then
            mkdir $archiveFolder
        fi

        currNum=0
        while [ -f $archiveFolder/oldDir.`basename $targetFolder`".$currNum.tgz" ]
        do
            currNum=$[$currNum+1]
        done
        tar -cvzf $archiveFolder/oldDir.`basename $targetFolder`".$currNum.tgz" $targetFolder
        rm -r $targetFolder
    fi
    mkdir $targetFolder
}


stringMergeDir=$stringtieDir"/stringMerged/"


if [ $runStringMerge == "yes" ]; then

    reuse UGER

    if [ ! -d $stringMergeDir ]; then
        mkdir $stringMergeDir
    fi

    for file in $stringtieDir/*/*_stringtie
    do
        echo $file >> $stringMergeDir/samples.manifest
    done

    qsub -cwd -b y -V -N stringMerge -l h_vmem=5G -j y -o $stringMergeDir/stringMerge.out $stringtieSoftwareDir/stringtie -x chrM --merge -G /seq/vgb/swofford/ref/Canis_familiaris.CanFam3.1.95.chrRenamed.gtf -o $stringMergeDir/stringtie_Merged.gtf $stringMergeDir/samples.manifest 

fi

abundanceDir="abundances/"
oldAbundanceDir="oldAbundances/"
cufflinksAbundanceDir=$abundanceDir/"cufflinks/"
stringtieAbundanceDir=$abundanceDir/"stringtie/"
stringQuantDir=$stringtieAbundanceDir"/initQuant"
quantDir=$cufflinksAbundanceDir"/initQuant/"


if [ $runStringQuant == "yes" ]; then

    mkdir $stringtieAbundanceDir
    mkdir $stringQuantDir
    for file in $alignmentDir/*/*bam
    do
	id=`basename $file | sed 's:.sorted.bam::'`
	mkdir $stringQuantDir/$id

	qsub -b y -N stringQuant_$id -o $stringQuantDir/$id/$id.out -cwd -j y -l h_vmem=10G -l h_rt=10:00:00 -V $stringtieSoftwareDir/stringtie -x chrM -e -B -G $stringMergeDir/stringtie_Merged.gtf -A $stringQuantDir/$id/"$id"_abund.tab -o $stringQuantDir/$id/$id.gtf $file
    done
fi

if [ $assemblyFASTA == "yes" ]; then
    reuse Cufflinks
    reuse UGER
    qsub -cwd -V -j y -o $stringMergeDir/ -b y gffread $stringMergeDir/stringtie_Merged.gtf -g $dogRef -w $stringMergeDir/stringtie_Merged.fasta
fi

if [ $runFastQC == "yes" ]; then
    cd ./alignments
    reuse .fastqc-0.11.4
    reuse Picard-Tools
    mkdir ./fastqc
    mkdir ./jobOut
    for file in */*sorted.bam
    do
	testing=$(readlink -m $file)
        qsub -b y -cwd -V -o ./jobOut/ fastqc --outdir ./fastqc/ $testing
    done
fi
