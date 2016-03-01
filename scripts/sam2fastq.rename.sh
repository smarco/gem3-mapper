#!/bin/bash

awk '
BEGIN{
  id=0;
  map["A"] = "T";
  map["C"] = "G";
  map["G"] = "C";
  map["T"] = "A";
  map["N"] = "N";
}
{
  if ($0 !~ /^@/) {
    if (and($2,0x10)) {
      printf("@%s.%07d\n",$1,id);
      for (i = length($10); i; i--) {
        printf "%s", map[substr($10, i, 1)]
      }
      printf("\n+\n");
      for (i = length($11); i; i--) {
        printf "%s", substr($11, i, 1)
      }
      printf("\n");
    } else {
      printf("@%s.%07d\n%s\n+\n%s\n",$1,id,$10,$11);
    }
    id++;
  }
}' $1
