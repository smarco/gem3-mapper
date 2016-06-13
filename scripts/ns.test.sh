#!/bin/bash

KEY=$1
ERROR=$2

./bin/gem-constructor -s edit-brute-supercondensed -i $KEY -n $ERROR> ns.brute
./bin/gem-constructor -s edit-partition -i $KEY -n $ERROR > ns.partition

DIFF=$(diff <(sort ns.brute | uniq ) <(sort ns.partition | uniq))
if [ "$DIFF" != "" ] 
then
  echo "NS FAILED"
  diff -y <(sort ns.brute | uniq ) <(sort ns.partition | uniq) | less
else
  echo "PASSED!" 
fi
