#!/bin/bash 

INPUT=$1

cat $INPUT | awk '{if ($0 !~ /^@/) print}' | awk '
BEGIN {
  tag="";
}
{
  if ($1!=tag) {
    tag=$1;
    print $1"\t"$3"\t"$4"\t"$5"\t"$6
  }
}' 
