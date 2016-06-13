#!/bin/bash

awk -v label=$1 -v nt=$2 '
BEGIN{
  num_reads=0;
  num_bases=0;
  limited=0;
}
{
  nr = (NR-1)%4;
  if (nr==0) {
    $1=sprintf("@%s.%010d",label,num_reads);
    print $1;
    num_reads++;
  } else if (nr==2) {
    print "+";
  } else {
    if (nr==1) {
      if (num_bases + length($0) > nt) {
        limited = nt - num_bases;
      }
      num_bases += length($0);
    }
    if (limited > 0) {
      a = substr($0,1,limited);
      print a;
    } else {
      print;
    }
    if (nr==3 && num_bases >= nt) {
      exit;
    }
  }
}'
