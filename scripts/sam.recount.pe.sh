#!/bin/bash

awk -v lab=$2 '
BEGIN {
  OFS="\t"
  id = 0;
  end = 1;
}
{
  $1=sprintf("Sim.Illumina.l%s.%010d",lab,id); 
  print;
  if (end==1) {
    end = 2;
  } else {
    end = 1;
    id++;
  } 
}' $1
