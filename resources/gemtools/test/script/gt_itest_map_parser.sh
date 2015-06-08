#!/bin/bash

BASH=/bin/bash
TIMEFORMAT=%R;
FILE_ERROR=gt_itest_map_parser.err
FILE_OUTPUT=gt_itest_map_parser.out
for file in ../datasets/*
do 
	# Log
	echo "#$file" >> $FILE_OUTPUT;
	echo "#$file" >> $FILE_ERROR;
	# Run tests
	FILE_SIZE=$(du -h $file | awk '{print $1}');
	echo "Parsing $file ($FILE_SIZE)"; 
    TIME_COUNT=$( { time wc $file > /dev/null; } 2>&1 );
    echo -n "    wc=$TIME_COUNT";
    \time $BASH -c "../bin/gt.stats -i $file -q 2>> $FILE_ERROR >> $FILE_OUTPUT" 2>&1 | \
      perl -e '($min,$sec,$micro)=<>=~/(\d+)\:(\d+)\.(\d+)elapsed/; ($mem)=<>=~/(\d+)minor/; $time=$min*60+$sec; print "    gt=$time.$micro ($mem pages)\n"';
done
echo ">>> Total Error: ";
# cat $FILE_ERROR;
echo ">>> Total Output: ";
# cat $FILE_OUTPUT;
rm $FILE_ERROR $FILE_OUTPUT;