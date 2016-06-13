#!/bin/bash

MASON=~smarco/bin/mason
REFERENCE=/scratch/077-hpca4se-bioinf/smarco/GEM3.GPU/references/hsapiens_v37.fa

NUM_READS=$1
READ_LENGTH=$2
LL=$3 # Default=375
LE=$4 # Default=100
PREFIX=$5

OUTPUT_FILE=Sim.Illumina.l${READ_LENGTH}.${PREFIX}
\time -v $MASON illumina --num-reads $NUM_READS -n $READ_LENGTH --output-file $OUTPUT_FILE.fastq -hn 2 -sq -ll $LL -le $LE --mate-pairs --read-naming 2 --very-verbose $REFERENCE > $OUTPUT_FILE.out 2> $OUTPUT_FILE.err

