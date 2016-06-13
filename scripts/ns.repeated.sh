#!/bin/bash

./bin/gem-constructor -s edit-partition -i $1 -n $2 2> tmp > /dev/null
sort tmp | awk '
BEGIN{
  tag="";
  rep="";
  c=0;
}
{
  if ($1==tag) {
    if (c==0) {
      print "---"
      print rep;
      print $0;
    } else {
      print $0;
    }
    rep=$0; 
    c++;
  } else {
    tag=$1;
    rep=$0;
    c=0;
  }
}' | less

rm tmp;
