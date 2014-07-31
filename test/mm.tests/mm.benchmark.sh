#!/bin/bash

cp ~/Indexes/human.gem $TMPDIR/human.gem

echo "==> gt_mm_bulk_malloc" >&2
./bin/gt.construct -s 0 -n 5070346304
./bin/gt.construct -s 0 -n 5070346304

echo "" >&2
echo "==> gt_mm_bulk_malloc_temp" >&2
./bin/gt.construct -s 1 -n 5070346304
./bin/gt.construct -s 1 -n 5070346304

echo "" >&2
echo "==> gt_mm_bulk_mmalloc PRIVATE" >&2
./bin/gt.construct -s 2 -n 5070346304 -1 0
./bin/gt.construct -s 2 -n 5070346304 -1 0

echo "" >&2
echo "==> gt_mm_bulk_mmalloc SHARED" >&2
./bin/gt.construct -s 2 -n 5070346304 -1 1
./bin/gt.construct -s 2 -n 5070346304 -1 1

flush.file.cache $TMPDIR/human.gem;
echo "" >&2
echo "==> gt_mm_bulk_mmap_file PRIVATE" >&2
./bin/gt.construct -i $TMPDIR/human.gem -s 3 -1 0; flush.file.cache $TMPDIR/human.gem;
./bin/gt.construct -i $TMPDIR/human.gem -s 3 -1 0; flush.file.cache $TMPDIR/human.gem;

echo "" >&2
echo "==> gt_mm_bulk_mmap_file SHARED" >&2
./bin/gt.construct -i $TMPDIR/human.gem -s 3 -1 1; flush.file.cache $TMPDIR/human.gem;
./bin/gt.construct -i $TMPDIR/human.gem -s 3 -1 1; flush.file.cache $TMPDIR/human.gem;

echo "-------------------------------------------------------------------------------" >&2
for i in 1 2 4 8 10 12
do
  echo "" >&2
  echo "==> gt_mm_bulk_load_file T=$i 256MB" >&2
  ./bin/gt.construct -i $TMPDIR/human.gem -s 4 -1 $i; flush.file.cache $TMPDIR/human.gem;
  ./bin/gt.construct -i $TMPDIR/human.gem -s 4 -1 $i; flush.file.cache $TMPDIR/human.gem;
done

echo "-------------------------------------------------------------------------------" >&2
for i in 1 2 4 8 10 12
do
  echo "" >&2
  echo "==> gt_mm_bulk_mload_file T=$i 256MB" >&2
  ./bin/gt.construct -i $TMPDIR/human.gem -s 5 -1 $i; flush.file.cache $TMPDIR/human.gem;
  ./bin/gt.construct -i $TMPDIR/human.gem -s 5 -1 $i; flush.file.cache $TMPDIR/human.gem;
done

echo "" >&2
echo "==> gt_mm_bulk_mmap_populate_file PRIVATE" >&2
./bin/gt.construct -i $TMPDIR/human.gem -s 6 -1 0; flush.file.cache $TMPDIR/human.gem;
./bin/gt.construct -i $TMPDIR/human.gem -s 6 -1 0; flush.file.cache $TMPDIR/human.gem;

echo "" >&2
echo "==> gt_mm_bulk_mmap_populate_file SHARED" >&2
./bin/gt.construct -i $TMPDIR/human.gem -s 6 -1 1; flush.file.cache $TMPDIR/human.gem;
./bin/gt.construct -i $TMPDIR/human.gem -s 6 -1 1; flush.file.cache $TMPDIR/human.gem;


