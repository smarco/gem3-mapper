#!/bin/bash

MODE=$1
KEY=$2
ERROR=$3

if [[ $MODE = "hamming" ]];
then
  ./bin/gem-constructor -s hamming-brute -i $KEY -n $ERROR> ns.brute
  ./bin/gem-constructor -s hamming-partition -i $KEY -n $ERROR > ns.partition
  # ./bin/gem-constructor -s hamming-regions -i $KEY -n $ERROR > ns.partition
else
  ./bin/gem-constructor -s edit-brute-supercondensed -i $KEY -n $ERROR> ns.brute
  # ./bin/gem-constructor -s edit-partition -i $KEY -n $ERROR > ns.partition
  ./bin/gem-constructor -s edit-regions -i $KEY -n $ERROR > ns.partition
fi

MISSING=$(comm -2 -3 <(sort ns.brute) <(sort ns.partition))
if [ "$MISSING" != "" ] 
then
  echo "NS FAILED"
  diff -y <(sort ns.brute | uniq ) <(sort ns.partition | uniq) | less
  ## comm -2 -3 <(sort ns.brute) <(sort ns.partition) | less
else
  COMMON=$(comm -1 -2 <(sort ns.brute) <(sort ns.partition) | wc -l)
  EXTRA=$(comm -1 -3 <(sort ns.brute) <(sort ns.partition) | wc -l)
  DUPS=$(sort ns.partition | uniq -d | wc -l)
  echo "PASSED! Total $COMMMON ($EXTRA extra) ($DUPS dups)" 
fi
