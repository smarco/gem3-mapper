#!/bin/bash

function refine () {
  awk '{nr=(NR-1)%4; if (nr==3) {printf $1"\n"} else if (nr==2) { printf "+\t"} else {printf $1"\t"}}' $1
}

paste <(refine $1) <(refine $2) | awk '
{
  print $1"\n"$2"\n"$3"\n"$4;
  print $5"\n"$6"\n"$7"\n"$8;
}'
