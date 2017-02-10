#!/bin/bash

VA=$1
VB=$2

for specie in chr1 human
do

for mapping in fast sensitive
do

for mode in se pe
do

echo -en "Checking sample.$specie.l100.$mode.$mapping\t"
DIFF=$(diff sample.$specie.l100.$mode.$mapping.$VA.map sample.$specie.l100.$mode.$mapping.$VB.map)
if [ "$DIFF" == "" ]; then echo "[OK]"; else echo "[FAILED]"; fi;

done
done
done
