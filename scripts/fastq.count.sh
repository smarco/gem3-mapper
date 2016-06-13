#!/bin/bash

awk '
BEGIN{nt=0;reads=0;}
{
  if ((NR-1)%4==1) {
    nt += length($1);
    reads++;
  }
}
END{print FILENAME": "reads" reads "nt" nt";}' $1
