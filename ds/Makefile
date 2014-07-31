SHELL=/bin/sh

CC=gcc

#these are for testing
#CFLAGS = -g -W -Wall -Winline -O2 

#these are for maximum speed
CFLAGS=-g -O3 -fomit-frame-pointer -W -Wall -Winline -m32 \
       -DDEBUG=0 -DNDEBUG=1  


.PHONY: all
all : ds unbwt bwt testlcp


# deep-shallow suffix sort algorithm
ds: suftest2.o ds_ssort.a 
	$(CC) $(CFLAGS) -o ds suftest2.o ds_ssort.a

# archive containing the ds sort algorithm
ds_ssort.a: globals.o ds.o shallow.o deep2.o helped.o blind2.o
	ar rc ds_ssort.a globals.o ds.o shallow.o deep2.o helped.o blind2.o

# archive containing the bwt and lcp auxiliary routines 
bwtlcp.a: bwt_aux.o lcp_aux.o
	ar rc bwtlcp.a bwt_aux.o lcp_aux.o

# compare several linear time lcp algorithms
testlcp: testlcp.c bwtlcp.a ds_ssort.a 
	 $(CC) $(CFLAGS) -o testlcp testlcp.c bwtlcp.a ds_ssort.a 

# inverse bwt
unbwt: unbwt.c
	 $(CC) $(CFLAGS) -o unbwt unbwt.c 

# bwt using ds_ssort
bwt: bwt.c ds_ssort.a
	 $(CC) $(CFLAGS) -o bwt bwt.c ds_ssort.a 

# pattern rule for all objects files
%.o: %.c *.h
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f *.o *.a






