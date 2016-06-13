#!/bin/bash

paste $1 $2 | awk '
function pdiff(p1,p2)
{
  if (p1 < p2) return (p2 - p1); else return (p1 - p2); 
}
{
  if ($1!=$6) {
    print "Wrong TAG Synch ("$0")" > "/dev/stderr";
    exit;
  }
  $6=""; 
  if ($2!=$7 || pdiff($3,$8) > 15) print
}'
