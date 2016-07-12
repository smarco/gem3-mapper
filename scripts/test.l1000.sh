#!/bin/bash

./bin/gem-mapper -I ../data/chr1.gem -i ../data/sample.chr1.l1000.lite.fastq -o sample.chr1.l1000.se.cpu.$1.map -F MAP
./bin/gem-mapper -I ../data/chr1.gem -i ../data/sample.chr1.l1000.lite.fastq -o sample.chr1.l1000.se.gpu.$1.map -F MAP --cuda=emulated
