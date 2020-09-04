#!/bin/bash
# Script to generate counts for CRISPR
# Requires following tools:
# - Bowtie
# Assumes each gRNA to be 20 bases

set -e
set -u
set -o pipefail

##check arguments
if [ "$#" -ne 10 ]; then
    printf "Usage: $(basename $0) [-1 full path to read1fastq (gz file only)] [-i Reference Fasta File] [-o Full path to Output Directory] [-3 bases to trim from 3' end] [-5 bases to trim from 5' end]\n";
    exit
fi

while getopts ":1:3:5:o:h:i:" OPTION; do
    case $OPTION in
        1)
            read1fastq=$OPTARG
            ;;
        3)
            trim3=$OPTARG
            ;;
        5)
            trim5=$OPTARG
            ;;
        i)
            fasta_file=$OPTARG
            ;;   
        o)
            outdir=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
	        echo "script usage: $(basename $0) [-1 full path to read1fastq (gz file only)] [-i Reference Fasta File] [-o Full path to Output Directory] [-3 bases to trim from 3' end] [-5 bases to trim from 5' end]\n" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG missing" >&2
	        echo "script usage: $(basename $0) [-1 full path to read1fastq (gz file only)] [-i Reference Fasta File] [-o Full path to Output Directory] [-3 bases to trim from 3' end] [-5 bases to trim from 5' end]\n" >&2
            exit 1
            ;;
        h)
            echo -e "-h print this screen and exit"
            exit $E_OPTERROR
            ;;
    esac
done
shift $((OPTIND -1))

#create output director if does not exit
if [ ! -d "$outdir" ]; then
  mkdir $outdir
else
  echo "Output Directory Already Exists"
  exit 2
fi

#get filename and extension
filename="${read1fastq##*/}"
fileext="${read1fastq#*.}"

##creates log file with input commands
echo "Log file $(date +%y_%m_%d)" > $outdir/log_file.txt;
echo "Command line input: "$0" "$@>>$outdir/log_file.txt;
echo "FASTQ file: $read1fastq" >>$outdir/log_file.txt;
echo "FASTA input: $fasta_file" >>$outdir/log_file.txt;
echo "Output file directory: $outdir" >>$outdir/log_file.txt;
echo "Input Filename: $filename" >>$outdir/log_file.txt;
echo "Extension: $fileext" >>$outdir/log_file.txt;

mkdir $outdir/bowtie_index
cp $fasta_file $outdir/bowtie_index/
cd $outdir/bowtie_index/

readLength=$(echo "$(zcat $read1fastq | head -n 4 | awk '{if(NR%4==2) print length($1)}')")

#bowtie commands, arguments
bowtie-build $(basename -- "$fasta_file") guides
cd ..

#read length
grnalength=`expr $readLength - $trim5 - $trim3`

if [ $grnalength -ne 20 ]; 
then
    echo "Read Length not consistent with having a 20 base gRNA sequence\n"
    exit 2
else
    #alignment
    bowtie -5 $trim5 -3 $trim3 -n 0 -l 20 -p 4 -S $outdir/bowtie_index/guides $read1fastq $outdir/${filename/.txt.gz/}_AlignedGuides.sam

    #counts
    grep -v ^\@ $outdir/${filename/.txt.gz/}_AlignedGuides.sam | awk -F"\t" '{ if($4=1 && $13=="MD:Z:20") print $3 }' | sort | uniq -c | awk '{ print $2"\t"$1 }' > $outdir/${filename/.txt.gz/}.count
fi
