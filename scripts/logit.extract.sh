#!/bin/sh

HEADER=$1
TP=$2
PREDICTORS=$3

# Compose
echo "ALL"
wc -l $TP $PREDICTORS
paste $TP $PREDICTORS > tmp
cat $HEADER tmp > logit.all
rm tmp

# Split
awk '{if (NR==1 || ($3=="signal")) print}' logit.all > logit.signal
awk '{if (NR==1 || ($3=="unique")) print}' logit.all > logit.uniq
awk '{if (NR==1 || ($3=="mmap")) print}' logit.all > logit.mmaps
awk '{if (NR==1 || ($3=="mmapD1")) print}' logit.all > logit.mmaps.d1
awk '{if (NR==1 || ($3=="tie")) print}' logit.all > logit.ties

echo "PARTS: logit.signal logit.uniq logit.mmaps logit.mmaps.d1 logit.ties"
wc -l logit.signal logit.uniq logit.mmaps logit.mmaps.d1 logit.ties
