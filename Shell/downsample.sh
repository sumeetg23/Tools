#!/bin/bash

# exit on errr
set -e

log() {
    >&2 echo $*
}

VERSION=1.0
log VERSION $VERSION

INDIR=$1
OUTDIR=$2
READS=$3
PATTERN=$4

FPATTERN=${PATTERN}'*_R1_*.fastq.gz'
RPATTERN=${PATTERN}'*_R2_*.fastq.gz'

[[ ! -e $OUTDIR ]] && mkdir $OUTDIR

# get first file name - forward
R1FILE=`ls -1 ${INDIR}/${FPATTERN} | head -1`
R2FILE=`ls -1 ${INDIR}/${RPATTERN} | head -1`

if [[ ! -e $R1FILE ]]; then
    error 'No forward strands found!'
    exit 1
fi
if [[ ! -e $R2FILE ]]; then
    error 'No reverse strands found!'
    exit 1
fi

R1FILE=`basename $R1FILE`
R2FILE=`basename $R2FILE`

# get a random number
SEED=$RANDOM

log 'Input files :'
ls ${INDIR}/${FPATTERN}
ls ${INDIR}/${RPATTERN}

SAMPLESIZE=$READS

log 'Running:'

# create a temp file to store the PID of seqtk
PIDFILE=`mktemp`

# create forward downsample file -- and put it in the background
COMMAND1="seqtk sample -s $SEED"
COMMAND2="gzip --to-stdout"
log "$COMMAND1 <(zcat ${INDIR}/${R1FILE}) $SAMPLESIZE | $COMMAND2 > ${OUTDIR}/${R1FILE} &"
( echo $BASHPID > $PIDFILE; exec $COMMAND1 <(zcat ${INDIR}/${R1FILE}) $SAMPLESIZE ) | $COMMAND2 > ${OUTDIR}/${R1FILE} &

# create reverse downsample file
log "$COMMAND1 <(zcat ${INDIR}/${R2FILE}) $SAMPLESIZE | $COMMAND2 > ${OUTDIR}/${R2FILE}"
$COMMAND1 <(zcat ${INDIR}/${R2FILE}) $SAMPLESIZE | $COMMAND2 > ${OUTDIR}/${R2FILE}