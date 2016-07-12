#!/bin/bash

SOURCE=$1
GEM=$2
BWA=$3

PSOURCE=prof.source
PGEM=prof.gem
PBWA=prof.bwa

# Make summaries
./scripts/sam.summary.se.sh $SOURCE > $PSOURCE.sum
./scripts/sam.summary.se.sh $GEM    > $PGEM.sum
./scripts/sam.summary.se.sh $BWA    > $PBWA.sum

# Make diff-s
./scripts/summary.diff.se.sh $PSOURCE.sum $PGEM.sum > $PGEM.diff
./scripts/summary.diff.se.sh $PSOURCE.sum $PBWA.sum > $PBWA.diff

# Make same-s
./scripts/summary.same.se.sh $PSOURCE.sum $PGEM.sum > $PGEM.same
./scripts/summary.same.se.sh $PSOURCE.sum $PBWA.sum > $PBWA.same

# Make lhits-s (left hits)
join $PGEM.diff $PBWA.same | awk '{$10="";$11="";$12="";$13="";print;}' | tr -s " " | sort -n -k12 > $PGEM.lhits
join $PBWA.diff $PGEM.same | awk '{$10="";$11="";$12="";$13="";print;}' | tr -s " " | sort -n -k12 > $PBWA.lhits

# Make lhits-s FASTQ
join <(sort $PGEM.lhits | awk '{print $1}') <(awk '{if ($0 !~ /^@/) print}' $GEM) | ./scripts/sam2fastq.se.sh > $PGEM.lhits.fastq

# Make
./scripts/summary.mapq.zero.sh $PGEM.sum $PBWA.sum > $PGEM.mapq.zero
