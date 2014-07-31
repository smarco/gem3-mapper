/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt.c 
   Ver 1.0    14-oct-02
   compute the bwt using the ds_ssort routines 
   The output file has the following format: 
      (bwt) N bytes + (eof position, LSB first) 4 bytes

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


/* --- read proptotypes and typedef for ds_ssort --- */
#include "ds_ssort.h"

// ----- global variables -------
static int Verbose;
FILE *Outfile;         // output file 
FILE *Infile;          // input file 
int Infile_size;       // size input file;


// ----- prototypes ----------
typedef unsigned char UChar;
void fatal_error(char *);
void out_of_mem(char *f);


/* ***************************************************************
   this procedures demostrates the use of the ds_ssort routine
   It computes the suffix array for the content of Infile and use
   it to compute the Burrows-Wheeler Transform, which is then
   written to Outfile.
   The original file can be retreived using the unbwt tool.
   *************************************************************** */
void bwt_file(void)
{
  UChar *text;
  int *sa, n, overshoot, eof_pos=0, i;

  // ----- init ds suffix sort routine
  overshoot=init_ds_ssort(500,2000);
  if(overshoot==0)
    fatal_error("ds initialization failed! (bwt_file)\n");
  // ----- allocate text and suffix array
  n = Infile_size;                               // length of input text
  sa=malloc((n)*sizeof *sa);                     // suffix array
  text=malloc((n+overshoot)*sizeof *text);       // text
  if (! sa || ! text) out_of_mem("bwt_file");
  // ----- read text and build suffix array
  rewind(Infile); 
  i=fread(text, (size_t) 1, (size_t) n, Infile);
  if(i!=n) fatal_error("Error reading the input file!");
  ds_ssort(text,sa,n);                           // sort suffixes
  // ----- write bwt to Outfile
  if(Verbose>1) fprintf(stderr,"Writing bwt to output file\n");
  putc(text[n-1],Outfile);              // L[0] = Text[n-1]
  for(i=0;i<n;i++) 
    if(sa[i]!=0)
      putc(text[sa[i]-1],Outfile);      // write char preceeding sa[i]
    else
      eof_pos = i;                      // store position of EOF
  for(i=0;i<4;i++) {
    putc(eof_pos & 0xFF, Outfile);      // write EOF position using 4 bytes
    eof_pos >>= 8;
  }
  free(text);
  free(sa);
}


/* ***************************************************************
   Tool for computing the bwt using the ds suffix sorter 
   *************************************************************** */
int main(int argc, char *argv[])
{  
  double getTime ( void );
  void open_files(char *infile_name, char *outfile_name);
  void bwt_file(void);
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
    fprintf(stderr,"{infile}.bwt\n\n");
    exit(0);
  }

  /* ------------- read command line options ----------- */
  Verbose=0;                       // be quiet by default 
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
    fprintf(stderr,"\n             bwt  Ver 1.0\n");
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
  bwt_file();
  end = getTime();
  if(Verbose)
    fprintf(stderr,"Elapsed time: %f seconds.\n", end-start);
  return 0;
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
    /* if outfile was not given add ".bwt" to infile_name */
    outfile_name = (char *) malloc(strlen(infile_name)+5);
    outfile_name = strcpy(outfile_name, infile_name);
    outfile_name = strcat(outfile_name, ".bwt");
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




