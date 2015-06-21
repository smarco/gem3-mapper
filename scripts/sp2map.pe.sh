#!/bin/bash
#
# CMD: ./sp2map rmap map (generates the corresponding MAP file)

RMAP=$1
MAP=$2

awk '{print $1}' $RMAP > t1

join t1 $MAP > ${RMAP%.rmap}.map

rm t1
