#!/bin/bash

VA=$1
VB=$2

for specie in chr1 human
do

for mapping in fast sensitive
do

for mode in se pe
do
  # CMP
  echo -en "Checking sample.$specie.l100.$mode.$mapping\t"
  DIFF=$(diff sample.$specie.l100.$mode.$mapping.$VA.map sample.$specie.l100.$mode.$mapping.$VB.map)
  if [ "$DIFF" == "" ]; then echo -n "[OK]"; else echo -n "[FAILED]"; fi;
  # TIME
  grep GEMMapper sample.$specie.l100.$mode.$mapping.$VA.log | awk '{printf " "$8}'
  grep GEMMapper sample.$specie.l100.$mode.$mapping.$VB.log | awk '{printf " vs "$8"\n"}'
done
done
done
