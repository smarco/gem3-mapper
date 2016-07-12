#!/bin/bash

RMAP=$1

join -1 1 -2 11 $RMAP logit.metrics | tr ' ' '\t' > ${RMAP%.rmap}.metrics
