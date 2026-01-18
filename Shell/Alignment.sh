#!/bin/bash
# Script to generate counts for RNA-Seq
# Requires following tools:
# - STAR
# - HTSeq Count
# - feature count
# - samtools
# Currently supports only single end rna-seq reads

# works on LSF cluster

set -e
set -u
set -o pipefail

##check arguments
if [ "$#" -ne 10 ]; then
    printf "Usage: $(basename $0) [-l read1fastq (gz file only)] [-i STAR Index] [-g GTF File] [-o Output Directory] [-s Stranded ("yes" or "no" or "reverse")]\n";
    exit
fi

while getopts ":1:s:g:o:h:i:" OPTION; do
    case $OPTION in
        1)
            read1fastq=$OPTARG
            ;;
        i)
            star_index=$OPTARG
            ;;
        g)
            gtf_file=$OPTARG #sorted
            ;;
   
        o)
            outdir=$OPTARG
            ;;

	s)
	    stranded=$OPTARG
	    ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
	    echo "script usage: $(basename $0) [-l read1fastq (gz file only)] [-s STAR Index] [-g GTF File] [-d Output Directory]" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG missing" >&2
	    echo "script usage: $(basename $0) [-l read1fastq (gz file only)] [-s STAR Index] [-g GTF File] [-d Output Directory]" >&2
            exit 1
            ;;
        h)
            echo -e "-h print this screen and exit"
            exit $E_OPTERROR
            ;;
    esac
done
shift $((OPTIND -1))

#check is strand is specified properly
if [ $stranded != "yes" ] && [ $stranded != "no" ] && [ $stranded != "reverse" ]; then
   echo "Strand not properly specified, Only 1 of the following options"
   echo "yes; no; reverse. 'reverse' means 'yes' with reversed strand"
   exit 1
fi

#Strand info for qualimap (qc report)
qualimapstr="strand-specific-reverse"
case $stranded in
   yes)
	   $qualimapstr = "strand-specific-forward"
	   ;;
   no)
	   $qualimapstr = "non-strand-specific"
	   ;;
esac

#create output director if does not exit
if [ ! -d "$outdir" ]; then
  mkdir $outdir
else
  echo "Output Directory Already Exists"
  exit 2
fi

#random number generator
RANDOM=$$
jobid=$RANDOM

#get filename and extension
filename="${read1fastq##*/}"
fileext="${read1fastq#*.}"

##creates log file with input commands
echo "Log file $(date +%y_%m_%d)" > $outdir/log_file.txt;
echo "Command line input: "$0" "$@>>$outdir/log_file.txt;
echo "FASTQ file: $read1fastq" >>$outdir/log_file.txt;
echo "STAR index: $star_index" >>$outdir/log_file.txt;
echo "GTF file: $gtf_file" >>$outdir/log_file.txt;
echo "Output file directory: $outdir" >>$outdir/log_file.txt;
echo "Input Filename: $filename" >>$outdir/log_file.txt;
echo "Extension: $fileext" >>$outdir/log_file.txt;

#bowtie commands, arguments
unqstr="${jobid}align"
bsub -o $outdir/AlignJob-%J.txt -J $unqstr "STAR --outSAMtype BAM SortedByCoordinate --genomeDir $star_index --runThreadN 4 --outFilterMultimapNmax 1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --readFilesIn $read1fastq --readFilesCommand zcat --outFileNamePrefix $outdir/${filename}_ > $outdir/${filename}.starout"

#index the bam file
idtrack=$unqstr
unqstr="${jobid}index"
bsub -o $outdir/IndexJob-%J.txt -w $idtrack -J $unqstr "samtools index $outdir/${filename}_Aligned.sortedByCoord.out.bam"

#computes expression
idtrack=$unqstr
unqstr="${jobid}fccount"
bsub -o $outdir/FCCountJob-%J.txt -w $idtrack -J $unqstr "featureCounts -s 0 -a $gtf_file -o $outdir/${filename}_featureCounts.count $outdir/${filename}_Aligned.sortedByCoord.out.bam"

unqstr="${jobid}htseqcount"
bsub -o $outdir/htseqCountJob-%J.txt -w $idtrack -J $unqstr "htseq-count -m union -f bam -s $stranded $outdir/${filename}_Aligned.sortedByCoord.out.bam $gtf_file  > $outdir/${filename}_HTSeqCounts.count"

#qc report
unqstr="${jobid}qualimap"
bsub -o $outdir/QualimapJob-%J.txt -w $idtrack -J $unqstr "qualimap rnaseq --java-mem-size=8G -bam $outdir/${filename}_Aligned.sortedByCoord.out.bam -gtf $gtf_file -outdir $outdir/${filename}_QualiMap2/ -oc $outdir/QualiMap2/${filename/.bam/_count} -outfile $outdir/${filename/.gzAligned.sortedByCoord.out.bam/}_ -outformat html --sequencing-protocol $qualimapstr"
