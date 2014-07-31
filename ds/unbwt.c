/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   unbwt.c 
   Ver 1.0    6-oct-02
   compute the inverse bwt. 
   The input file must be in the following format: 
      (bwt) N bytes + (eof position, LSB first) 4 bytes
   The procedure currently uses 5N bytes.    

   Copyright (C) 2002 Giovanni Manzini (manzini@mfn.unipmn.it)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

   See COPYRIGHT file for further copyright information	   
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <assert.h>


typedef unsigned char UChar;
// ---- struct containing the (uncompressed) bwt 
typedef struct {
  UChar *bwt;
  int size;
  int eof_pos;
} bwt_data;


// ----- global variables -------
int Verbose=0;         // default: quiet mode
FILE *Outfile;         // output file 
FILE *Infile;          // input file 
int Infile_size;       // size input file;


// ----- prototypes ----------
void fatal_error(char *);
void out_of_mem(char *f);


/* ***************************************************************
   Tool for inverting the bwt computed by the ds suffix sorter 
   *************************************************************** */
int main(int argc, char *argv[])
{  
  double getTime ( void );
  void open_files(char *infile_name, char *outfile_name);
  void unbwt_file(void);
  double end, start;
  char *infile_name, *outfile_name;
  int c;
  extern char *optarg;
  extern int optind, opterr, optopt;

  if(argc<2) {
    fprintf(stderr, "Usage:\n\t%s infile ",argv[0]);
    fprintf(stderr,"[-v][-o outfile] \n");
    fprintf(stderr,"\t-v           verbose output\n");
    fprintf(stderr,"\t-o outfile   output file\n");
    fprintf(stderr,"If no output file name is given default is ");
    fprintf(stderr,"{infile}.y\n\n");
    exit(0);
  }

  /* ------------- read command line options ----------- */
  Verbose=0;
  infile_name=outfile_name=NULL;
  opterr = 0;
  while ((c=getopt(argc, argv, "do:n:v")) != -1) {
    switch (c)
    {
      case 'o':
        outfile_name = optarg; break;
      case 'v':
	Verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
    }
  }
  if(optind<argc)
     infile_name=argv[optind];

  // ----- check input parameters--------------
  if(infile_name==NULL) fatal_error("The input file name is required\n");

  if(Verbose>2) {
    fprintf(stderr,"\n*****************************************************");
    fprintf(stderr,"\n             unbwt  Ver 1.0\n");
    fprintf(stderr,"Created on %s at %s from %s\n",__DATE__,__TIME__,__FILE__);
    fprintf(stderr,"*****************************************************\n");
  }
  if(Verbose>1) {
    fprintf(stderr,"Command line: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n");
  }

  // ---- do the work ---------------------------
  start = getTime();
  open_files(infile_name,outfile_name);
  unbwt_file();
  end = getTime();
  if(Verbose)
    fprintf(stderr,"Elapsed time: %f seconds.\n", end-start);
  return 0;
}


/* ***********************************************************
   read the bwt from Infile, invert it and write the result to Outfile 
   *********************************************************** */
void unbwt_file(void)
{
  void read_bwt(FILE *f, bwt_data *b);
  void compute_rn_map(bwt_data *b, int *rn);
  void write_inverse_bwt_rn(bwt_data *b, int *rn, FILE *f);
  bwt_data bb, *b=&bb;
  int *rank_next;

  read_bwt(Infile,b);
  if(Verbose>1) fprintf(stderr,"bwt successfully read!\n");
  rank_next = (int *) malloc(b->size*sizeof(int *));
  if(rank_next==NULL)  out_of_mem("unbwt_file");
  compute_rn_map(b, rank_next);
  if(Verbose>1) fprintf(stderr,"rank_next map computed!\n");
  write_inverse_bwt_rn(b,rank_next,Outfile);
  free(rank_next);
  free(b->bwt);
}


// use bwt and rank_next map to compute the inverse
// bwt and write it to f 
void write_inverse_bwt_rn(bwt_data *b, int *rank_next, FILE *f)
{
  int i,j, written=0;

  // in the following, i is an index in the last column (the bwt)
  // and j is an index in the first column
  j = b->eof_pos;   // this is the position of the first text char 
  while(1) {      
    // get last column index corresponding to j
    i = rank_next[j];
    assert(i>=0 && i<b->size);
    putc(b->bwt[i],f);  // output char
    if(++written==b->size) break;  // if all chars have been written stop
    // now we move to the first column staying in the same row
    // hence we get j such that f[j] is the character following l[i]
    if(i<=b->eof_pos) j=i-1;
    else j=i;
  }
  // check that the last written char was bwt[0]
  if(i!=0) 
    fatal_error("Error writing inverse bwt (write_inverse_bwt_rn)\n");  
}


// read bwt from file f and store it into b
void read_bwt(FILE *f, bwt_data *b)
{
  int i,c,fsize;

  fseek(f,0,SEEK_END);
  fsize = ftell(f);
  if(fsize<5) 
    fatal_error("Invalid bwt file (read_bwt)\n");
  // compute size and alloc bwt array
  b->size = fsize-4;
  if((b->bwt=(UChar *)malloc(b->size))==NULL)
    out_of_mem("read_bwt");
  // read bwt
  rewind(f); 
  i=fread(b->bwt, (size_t) 1, (size_t) b->size, f);
  if(i!=b->size) 
    fatal_error("Error reading bwt file (read_bwt)\n");
  // read eof pos
  b->eof_pos=0;
  for(i=0;i<4;i++) {
    c=fgetc(f);
    if(c==EOF) 
      fatal_error("Error reading bwt file (read_bwt)\n");
    b->eof_pos |= (c<<(8*i));
  }
  if(b->eof_pos<0 || b->eof_pos>=b->size) 
    fatal_error("Illegal eof in bwt file (read_bwt)\n");
}


// compute the rank_next map  and store it into rn[] 
// which must be an already allocated array of size b->size 
// rn[i] is the position in the bwt (the last column of the BWT matrix)
// of the character f[i] (the first column of the BWT matrix)
// ignoring the eof symbol
void compute_rn_map(bwt_data *b, int *rn)
{
  int i, occ[256];

  // clear occ
  for(i=0;i<256;i++) occ[i]=0;
  // compute occ
  for(i=0;i<b->size;i++)
    occ[b->bwt[i]]++;
  // accumulate occ
  for(i=1;i<256;i++)
    occ[i] += occ[i-1];
  // shift occ
  assert(occ[255]==b->size);
  for(i=255;i>0;i--)
    occ[i] = occ[i-1];
  occ[0]=0;
  /* ----- compute rank_next mapping -------- */
  assert(rn!=NULL);
  for(i=0;i<b->size;i++)
    rn[occ[b->bwt[i]]++] = i;
}


/* ***************************************************************
   open input and output files; initializes Outfile, Infile, Infile_size
   *************************************************************** */
void open_files(char *infile_name, char *outfile_name)
{
  FILE *fopen(const char *path, const char *mode);

  /* ------ open input and output files ------ */
  if(infile_name==NULL)
    fatal_error("Pleaase provide the input file name (open_files)\n");
  else {
    Infile=fopen( infile_name, "rb"); // b is for binary: required by DOS
    if(Infile==NULL) {
      perror(infile_name);
      exit(1);
    }
  }
  // ---- store input file length in Infile_size
  if(fseek(Infile,0,SEEK_END)!=0) {
    perror(infile_name); exit(1);
  }
  Infile_size=ftell(Infile);
  if (Infile_size==-1) {
    perror(infile_name); exit(1);
  } 
  if (Infile_size==0) fatal_error("Input file empty (open_files)\n");

  // --- open output file
  if(outfile_name==NULL) {
    /* if outfile was not given add ".y" to infile_name */
    outfile_name = (char *) malloc(strlen(infile_name)+3);
    outfile_name = strcpy(outfile_name, infile_name);
    outfile_name = strcat(outfile_name, ".y");
  }
  Outfile = fopen( outfile_name, "wb"); // b is for binary: required by DOS
  if(Outfile==NULL) {
    perror(outfile_name);
    exit(1);
  }
}


void fatal_error(char *s)
{
  fprintf(stderr,"%s",s);
  exit(1);
}

void out_of_mem(char *f)
{
  fprintf(stderr, "Out of memory in function %s!\n", f);
  exit(1);
}

double getTime ( void )
{
   double usertime,systime;
   struct rusage usage;

   getrusage ( RUSAGE_SELF, &usage );

   usertime = (double)usage.ru_utime.tv_sec +
     (double)usage.ru_utime.tv_usec / 1000000.0;

   systime = (double)usage.ru_stime.tv_sec +
     (double)usage.ru_stime.tv_usec / 1000000.0;

   return(usertime+systime);
}







