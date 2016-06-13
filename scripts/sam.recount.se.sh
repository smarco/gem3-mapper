#!/bin/bash

awk -v lab=$2 '
BEGIN {
  OFS="\t"
  id = 0;
}
{
  $1=sprintf("Sim.Illumina.l%s.%010d",lab,id); 
  print;
  id++;
}' $1
