#!/bin/bash

awk '
BEGIN {
  OFS="\t"
  id = 0;
  end = 1;
}
{
  $1=sprintf("Sim.Illumina.l100.%010d",id); 
  print;
  if (end==1) {
    end = 2;
  } else {
    end = 1;
    id++;
  } 
}' $1
