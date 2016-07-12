#!/bin/bash

awk '
BEGIN{nr=0}
{
  if ($0 !~ /^@/) {
    if (nr%2==0) {
      printf $1"\t"
    } else {
      printf "\t"
    }
    printf $3"\t"$4"\t"$5;
    if (nr%2==1) printf "\n";
    nr++;
  }
}' $1
