/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   testlcp.c 
   Ver 1.0    9-dec-03
   Ver 1.1   28-apr-04     (uses bwt_aux)
   Ver 1.2   28-may-04     (uses lcp_aux) 
   Test algorithms for computing the lcp array. 

   Copyright (C) 2004 Giovanni Manzini (manzini@mfn.unipmn.it)

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
#include "bwt_aux.h"
#include "lcp_aux.h"

#define ALPHA_SIZE _BW_ALPHA_SIZE
typedef unsigned long long uint64;

// ----- global variables -------
static int Check_lcp;
static int Verbose;
static FILE *Infile;            // input file 
static int Infile_size;         // size input file;

// ----- prototypes ----------
void fatal_error(char *);
void out_of_mem(char *f);
double getTime ( void );

/* ***************************************************************
   Sample computation of the lcp array using the deep/shallow
   suffix sorter and the routines in the bwtlcp.a
   
   This procedure reads the input parameters from the command line
   and calls lcp_file() that does all the testing.
   *************************************************************** */
int main(int argc, char *argv[])
{  
  double getTime ( void );
  void open_files(char *infile_name);
  void lcp_file(void);
  double end, start;
  char *infile_name;
  int c;
  extern char *optarg;
  extern int optind, opterr, optopt;

  if(argc<2) {
    fprintf(stderr, "Usage:\n\t%s infile ",argv[0]);
    fprintf(stderr,"[-Cv]\n");
    fprintf(stderr,"\t-C           check lcp array (very slow!)\n");
    fprintf(stderr,"\t-v           verbose output\n\n");
    exit(0);
  }

  /* ------------- read command line options ----------- */
  Verbose=0;                       // be quiet by default
  Check_lcp=0;                     // do not check lcp by default
  infile_name=NULL;
  opterr = 0;
  while ((c=getopt(argc, argv, "Cv")) != -1) {
    switch (c)
    {
      case 'v':
	Verbose++; break;
      case 'C':
	Check_lcp=1; break;
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
    fprintf(stderr,"\n             testlcp  Ver 1.2\n");
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
  open_files(infile_name);
  lcp_file();
  end = getTime();
  if(Verbose)
    fprintf(stderr,"Elapsed time: %f seconds.\n", end-start);
  return 0;
}



/* ***************************************************************
   compute the lcp array using different linear time algorithms
   *************************************************************** */
void lcp_file(void)
{
  void cmp_lcp_array(int *lcp1, char *name1, int *lcp2, char *name2, int n);
  void check_lcp_array(uint8 *t, int n, int *sa, int *lcp);
  uint8 *text;
  int *sa, n, overshoot, i, extra_bytes;
  int occ[ALPHA_SIZE], *lcp9, *lcp=NULL;
  bwt_data b;
  double start, end;
  uint64 tot_lcp;

  // ----- init ds suffix sort routine
  overshoot=init_ds_ssort(500,2000);
  if(overshoot==0)
    fatal_error("ds initialization failed! (lcp_file)\n");
  // ----- allocate text and suffix array
  n = Infile_size;                               // length of input text
  sa=malloc((n+1)*sizeof *sa);                   // suffix array
  text=malloc((n+overshoot)*sizeof *text);       // text
  if (! sa || ! text) out_of_mem("lcp_file");
  // ----- read text and build suffix array
  rewind(Infile); 
  i=fread(text, (size_t) 1, (size_t) n, Infile);
  if(i!=n) fatal_error("Error reading the input file!");
  fprintf(stdout,"File size: %d\n",n);
  // ----- build suffix array ----------------
  start = getTime();
  ds_ssort(text,sa+1,n);                         // sort suffixes
  end=getTime();
  fprintf(stdout,"Suffix array construction: %.2f\n",end-start);

  // ---- compute lcp using 9n algorithm  
  start = getTime();
  for(i=0;i<ALPHA_SIZE;i++) occ[i]=0;
  for(i=0;i<n;i++) occ[text[i]]++;
  lcp9 = _lcp_sa2lcp_9n(text,n,sa,occ);
  end=getTime();
  fprintf(stdout,"lcp9 construction: %.2f\n",end-start);
  // check (very slowly) the correctness of lcp array
  if(Check_lcp) check_lcp_array(text,n,sa,lcp9);

  // --- compute lcp using 13n algorithm 
  start = getTime();
  lcp = _lcp_sa2lcp_13n(text, n, sa);
  end=getTime();
  fprintf(stdout,"lcp13 construction: %.2f\n",end-start);
  // ---- check lcp13 vs lcp9
  cmp_lcp_array(lcp, "lcp13", lcp9, "lcp9", n);
  free(lcp);

  // ---- compute lcp using 9.125 n algorithm (Makinen)  
  start = getTime();
  for(i=0;i<ALPHA_SIZE;i++) occ[i]=0;
  for(i=0;i<n;i++) occ[text[i]]++;
  lcp = _lcp_sa2lcp_9125n(text,n,sa);
  end=getTime();
  fprintf(stdout,"lcp9125 construction: %.2f\n",end-start);
  // ---- check lcp9125 vs lcp9
  cmp_lcp_array(lcp, "lcp9125", lcp9, "lcp9", n);
  free(lcp);


  // ---- compute lcp using 6n algorithm
  start = getTime();
  if( (b.bwt = (uint8 *) malloc(n+1)) == NULL)
    out_of_mem("lcp_file");
  _bw_sa2bwt(text, n, sa, &b);
  for(i=0;i<ALPHA_SIZE;i++) occ[i]=0;
  for(i=0;i<n;i++) occ[text[i]]++;
  extra_bytes = _lcp_sa2lcp_6n(text,&b,sa,occ);
  end=getTime();
  fprintf(stdout,"lcp6 construction: %.2f\n",end-start);
  fprintf(stdout,"Total memory for lcp6: %.2fn bytes\n",
	  6+(4.0*extra_bytes)/n);
  // ---- check lcp6 vs lcp9
  cmp_lcp_array(sa, "lcp6", lcp9, "lcp9", n);

  // ---- compute average lcp (which is now in sa[])
  tot_lcp=0;
  for(i=1;i<n;i++)
    tot_lcp += sa[i];
  fprintf(stdout,"Average lcp: %.2f\n",((double) tot_lcp)/(n-1));
  free(lcp9);
  free(text);
  free(sa);
}



/* *******************************************************************
   compare the lcp arrays cmputed by algs. name1, name2
   ******************************************************************* */ 
void cmp_lcp_array(int *lcp1, char *name1, int *lcp2, char *name2, int n)
{
  int i, diff;

  diff=0;
  for(i=2;i<=n;i++)
    if(lcp1[i]!=lcp2[i]) { 
      diff++; 
      if(Verbose>1)
	fprintf(stdout,"%s[%d]=%d  %s[%d]=%d\n",
		name1,i,lcp1[i],name2,i,lcp2[i]);
    }
  if(diff>0)
    fprintf(stdout,"FATAL ERROR! %s/%s differences: %d\n",name1,name2,diff);
}

/* *******************************************************************
   check the correctness of a lcp array with a (very slow) 
   character by character comparison of consecutive suffixes
   ******************************************************************* */ 
void check_lcp_array(uint8 *t, int n, int *sa, int *lcp)
{
  int i,j,k,h;

  for(i=2;i<n;i++) {
    j = sa[i-1]; k=sa[i];
    for(h=0;j+h<n && k+h<n;h++)
      if(t[j+h]!=t[k+h]) break;
    if(lcp[i]!=h) {
      fprintf(stdout,"FATAL ERROR! Incorrect LCP value: lcp[%d]=%d!=%d\n",
	      i,lcp[i],h);
      return;
    }
  }
  fprintf(stdout,"LCP array OK!\n");
}
 

/* ***************************************************************
   open input and output files; initializes Outfile, Infile, Infile_size
   *************************************************************** */
void open_files(char *infile_name)
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

