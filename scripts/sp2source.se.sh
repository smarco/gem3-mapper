#!/bin/bash

RMAP=$1

join $1 ../data/sample.chr1.l100.source.se.sam | tr ' ' '\t' > ${RMAP%.rmap}.source
