#!/bin/sh

echo "ALL"
wc -l logit.tp logit.metrics
paste logit.tp logit.metrics > tmp
cat logit.header tmp > logit.all
rm tmp

# awk '{if (NR==1 || ($10+$11>1 && ($2==$3 || $4==$5 || $6==$7) )) print }' logit.all > logit.ties
# awk '{if (NR==1 || ($10+$11>1 && !($2==$3 || $4==$5 || $6==$7) )) print}' logit.all > logit.mmaps
# awk '{if (NR==1 || ($10==1 && $11==0)) print }' logit.all > logit.uniq

# echo "MAIN: logit.all"
# wc -l logit.all
# echo "PARTS: logit.uniq logit.mmaps logit.ties"
# wc -l logit.uniq logit.mmaps logit.ties

awk '{if (NR==1 || (199<=$14 && $14<=250)) print}' logit.all > logit.uniq
awk '{if (NR==1 || (139<=$14 && $14<=190)) print}' logit.all > logit.mmaps
awk '{if (NR==1 || (9<=$14 && $14<=60)) print}' logit.all > logit.ties

echo "PARTS: logit.uniq logit.mmaps logit.ties"
wc -l logit.uniq logit.mmaps logit.ties
