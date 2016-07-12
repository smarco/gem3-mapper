#!/bin/bash

FILE=$1

awk -v f1=${FILE%.fastq}.1.fastq -v f2=${FILE%.fastq}.2.fastq '{
  nr=(NR-1)%8;
  if (nr <= 3) {
    print > f1
  } else {
    print > f2
  }
}' $FILE
