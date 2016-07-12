#!/bin/bash

awk '
BEGIN{
  id=0;
  end=1;
  map["A"] = "T";
  map["C"] = "G";
  map["G"] = "C";
  map["T"] = "A";
  map["N"] = "N";
}
{
  if ($0 !~ /^@/) {
    if (and($2,0x10)) {
      printf("@Sim.Illumina.l100.%010d/%d\n",id,end);
      for (i = length($10); i; i--) {
        printf "%s", map[substr($10, i, 1)]
      }
      printf("\n+\n");
      for (i = length($11); i; i--) {
        printf "%s", substr($11, i, 1)
      }
      printf("\n");
    } else {
      printf("@Sim.Illumina.l100.%010d/%d\n%s\n+\n%s\n",id,end,$10,$11);
    }
    if (end==1) {
      end=2;
    } else {
      end=1;
      id++;
    }
  }
}' $1
