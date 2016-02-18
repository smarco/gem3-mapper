#!/bin/bash 

INPUT=$1

cat $INPUT | awk '{if ($0 !~ /^@/) print}' | awk '
BEGIN {
  tag="";
  count=0;
}
{
  if (tag==$1) {
    if (count==1) {
      p2=$4
      map2=$3"\t"$4"\t"$5"\t"$6;
      count++;
    } 
    if (count==2) {
      if (paired==1) { printf "P\t" } else { printf "U\t" }
      if (p1 < p2) {
        print map1"\t"map2;
      } else {
        print map2"\t"map1;
      }
      count++;
    }
  } else {
    # Keep new TAG
    tag=$1;
    count=1;
    # Store 
    if (and($2,0x2)) { paired=1 } else { paired=0 }
    p1=$4
    map1=$3"\t"$4"\t"$5"\t"$6;
    # Print read-name
    printf $1"\t"
  }
}' 
