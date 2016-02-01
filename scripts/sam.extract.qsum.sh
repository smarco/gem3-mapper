#!/bin/bash

## NO cat $1 | awk '{if (and($2,0x2)) print}' | awk '{nr=(NR+1)%2; if(nr==1) printf $1" "$5; else printf " "$5"\n"}' > gem.pe.sum

cat ../data/sample.chr1.l100.bwa.pe.sam | awk '{if ($0 !~ /^@/) print}' | awk 'BEGIN{tag="";}{if (tag==$1) { printf " "$5 } else { printf "\n"$1" "$2" "$5;} tag=$1;}' | awk '{if (and($2,0x2)) print}'

awk '{if ($2==0) print}' gem.pe.sum > gem.pe.sum.q0
#awk '{if ($2==1) print}' gem.pe.sum > gem.pe.sum.q1
#awk '{if ($2==2) print}' gem.pe.sum > gem.pe.sum.q2
#awk '{if ($2==3) print}' gem.pe.sum > gem.pe.sum.q3

join gem.pe.sum.q0 bwa.pe.sum > gem.pe.sum.join.q0
#join gem.pe.sum.q1 bwa.pe.sum > gem.pe.sum.join.q1
#join gem.pe.sum.q2 bwa.pe.sum > gem.pe.sum.join.q2
#join gem.pe.sum.q3 bwa.pe.sum > gem.pe.sum.join.q3

awk '{print $4"-"$5}' join.q0 | sort | uniq -c | sort -n > gem.pe.sum.join.q0.sum
#awk '{print $4"-"$5}' join.q1 | sort | uniq -c | sort -n > gem.pe.sum.join.q1.sum
#awk '{print $4"-"$5}' join.q2 | sort | uniq -c | sort -n > gem.pe.sum.join.q2.sum
#awk '{print $4"-"$5}' join.q3 | sort | uniq -c | sort -n > gem.pe.sum.join.q3.sum
