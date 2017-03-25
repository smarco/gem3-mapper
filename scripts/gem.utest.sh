#!/bin/bash

# Arguments:: USAGE
if [[ "$#" -eq 0 ]]; then
  echo "USAGE:: ./gem.test.bash v0 se.fast/all human/chr1 cpu/gpu";
  exit;
else
  VERSION=$1;
fi
# Arguments:: MODE
if [[ "$#" -eq 1 ]]; then
  echo "Performing ALL ($VERSION)";
  MODE="all";
else
  echo "Performing $2 ($VERSION)";
  MODE=$2;
fi
# Arguments:: ORGANISM
if [[ "$#" -eq 2 ]]; then
  echo "Using ORGANISM:Chr1";
  ORGANISM="chr1";
  INDEX="../data/chr1.gem"
else
  if [[ $3 == "human" ]]; then
    echo "Using ORGANISM:Human";
    ORGANISM="human";
    INDEX="../data/hsapiens_v37.s32.gem"
  else
    echo "Using ORGANISM:Chr1";
    ORGANISM="chr1";
    INDEX="../data/chr1.gem"
  fi
fi
# Arguments:: PLATFORM
if [[ "$#" -eq 3 ]]; then
  echo "Using PLATFORM:CPU"
  PLATFORM=""
else
  if [[ $4 == "gpu" ]]; then
    echo "Using PLATFORM:GPU"
    PLATFORM="--gpu"
  else
    echo "Using PLATFORM:CPU"
    PLATFORM=""
  fi
fi

# Parameters
MAPPER=./bin/gem-mapper
MAPSET=./resources/gemtools/bin/gt.mapset
INPUT_SE=../data/sample.$ORGANISM.l100.se.fastq
INPUT_PE=../data/sample.$ORGANISM.l100.pe.fastq
SOURCE_SE=../data/sample.$ORGANISM.l100.source.se.sam
SOURCE_PE=../data/sample.$ORGANISM.l100.source.pe.sam

# SE fast
if [[ $MODE == "se.fast" || $MODE == "all" ]]; then
echo "GEMMapping SE Fast..."
PREFIX=sample.$ORGANISM.l100.se.fast.$VERSION
\time -v $MAPPER -I $INDEX -i $INPUT_SE -o $PREFIX.map -F MAP --mapping-mode=fast $PLATFORM 2> $PREFIX.log
$MAPSET -O specificity-profile --i1 $SOURCE_SE --i2 $PREFIX.map -o $PREFIX 2> $PREFIX.sp
fi

# SE sensitive
if [[ $MODE == "se.sensitive" || $MODE == "all" ]]; then
echo "GEMMapping SE Sensitive..."
PREFIX=sample.$ORGANISM.l100.se.sensitive.$VERSION
\time -v $MAPPER -I $INDEX -i $INPUT_SE -o $PREFIX.map -F MAP --mapping-mode=sensitive $PLATFORM 2> $PREFIX.log
$MAPSET -O specificity-profile --i1 $SOURCE_SE --i2 $PREFIX.map -o $PREFIX 2> $PREFIX.sp
fi

# PE fast
if [[ $MODE == "pe.fast" || $MODE == "all" ]]; then
echo "GEMMapping PE Fast..."
PREFIX=sample.$ORGANISM.l100.pe.fast.$VERSION
\time -v $MAPPER -I $INDEX -i $INPUT_PE -p -o $PREFIX.map -F MAP --mapping-mode=fast $PLATFORM 2> $PREFIX.log
$MAPSET -O specificity-profile --i1 $SOURCE_PE --i2 $PREFIX.map -p -o $PREFIX 2> $PREFIX.sp
fi

# PE sensitive
if [[ $MODE == "pe.sensitive" || $MODE == "all" ]]; then
echo "GEMMapping PE Sensitive..."
PREFIX=sample.$ORGANISM.l100.pe.sensitive.$VERSION
\time -v $MAPPER -I $INDEX -i $INPUT_PE -p -o $PREFIX.map -F MAP --mapping-mode=sensitive $PLATFORM 2> $PREFIX.log
$MAPSET -O specificity-profile --i1 $SOURCE_PE --i2 $PREFIX.map -p -o $PREFIX 2> $PREFIX.sp
fi

rm sp.*

