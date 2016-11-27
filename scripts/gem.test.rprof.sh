#!/bin/bash

INDEX=../data/chr1.gem
INPUT=../data/sample.chr1.l100.se.fastq
SOURCE=../data/sample.chr1.l100.source.se.sam

THRESHOLD=20
MAX_STEPS=4
ALPHABET=2

# Setup Prefix
OUTPUT=test.$THRESHOLD.$MAX_STEPS.$ALPHABET
# Run mapper
./bin/gem-mapper -I $INDEX -i $INPUT -o $OUTPUT.map -F MAP -t 1 --region-model-lightweight=$THRESHOLD,$MAX_STEPS,$ALPHABET,2 --profile 2> $OUTPUT.profile
# Run Specificity Profile
./resources/gemtools/bin/gt.mapset -C specificity-profile --i1 $SOURCE --i2 $OUTPUT.map -o $OUTPUT 2> $OUTPUT.sp
# Extract Summary
echo $OUTPUT" " > $OUTPUT.sum
grep -m1 "Num.Regions" $OUTPUT.profile | awk '{printf "("$3","}' >> $OUTPUT.sum
grep -m1 "Candidate.Positions" $OUTPUT.profile | awk '{printf $3","}' >> $OUTPUT.sum
grep -m1 "Candidate.Regions" $OUTPUT.profile | awk '{printf $3";"}' >> $OUTPUT.sum
grep -m1 "RANKS.Mapper" $OUTPUT.profile | awk '{printf $3")"}' >> $OUTPUT.sum
grep -m1 "Num.Mapped" $OUTPUT.sp | awk '{printf "("$4","}' >> $OUTPUT.sum
grep -m1 "Num.Reads.Any.TP" $OUTPUT.sp | awk '{printf $4","}' >> $OUTPUT.sum
grep -m1 "Num.Reads.First.TP" $OUTPUT.sp | awk '{printf $4")"}' >> $OUTPUT.sum

awk '{if ($0 ~ /^60,/) printf " "$1 }' >> $OUTPUT.sum
awk '{if ($0 ~ /^40,/) printf " "$1 }' >> $OUTPUT.sum
awk '{if ($0 ~ /^20,/) print " "$1 }' >> $OUTPUT.sum
# Delete MAP
rm $OUTPUT.map
