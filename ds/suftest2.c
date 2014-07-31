/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   suftest2.c
   Written by G. Manzini using the suftest.c code 
   by N. Jesper Larsson (Copyright N. Jesper Larsson 1999) 
   
   Program to test a suffix sorting function. Reads a sequence of bytes from
   a file and calls the deep/shallow suffix sorting algorithm.
   For further information on deep/shallow suffix sorting look at the paper:
     G. Manzini, P. Ferragina
     Engineering a Lightweight Suffix Array Construction Algorithm.
     Proc. 10th European Symposium on Algorithms (ESA '02).
     Springer Verlag Lecture Notes in Computer Science n. 2461, pp 698-710.
     http://www.mfn.unipmn.it/~manzini/lightweight

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

   See README file for further copyright information
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <sys/times.h>
#include "common.h"


/* ---- global variables (modifiable by command line options) ----- */
extern int Anchor_dist;                // distance between anchors
extern int Shallow_limit;              // limit for shallow_sort
extern int _ds_Verbose;                // how verbose it the algorithm?
extern int _ds_Word_size;              // # of bytes in word in mkqs
extern int Mk_qs_thresh;               // recursion limit for mk quicksort:
extern Int32 Max_pseudo_anchor_offset; // maximum offset considered when 
                                       // searching a pseudo anchor
extern Int32 B2g_ratio;                // maximum ratio bucket_size/group_size
                                       // accepted for pseudo anchor_sorting
extern Int32 Update_anchor_ranks;      // if!=0 update anchor ranks when 
                                       // determining rank for pseudo-sorting
extern Int32 Blind_sort_ratio;         // blind sort is used for groups of 
                                       // size <= Text_size/Blind_sort_ratio



#define MAX_LCP_SIZE 1000000    // maximum size of a recordable LCP


/* --------- prototypes ----------------- */ 
void ds_ssort(UChar *x, int *p, int n);
clock_t times(struct tms *buffer);
int check_global_variables(void);
void set_global_variables(void);
int compute_overshoot(void);


int main(int argc, char *argv[])
{
  void write_sa(char *filename, int *p, int n);
  void write_lcp(char *filename, UChar *x,int *p, int n);
  void write_bwt(char *filename, UChar *x,int *p, int n);
  void check_sa_ordering(UChar *x,int *p, int n, int);
  void print_sa_onscreen(UChar *x,int *p, int n, int);
  int c, *p, n;
  int print_sa, check_sa, num_opt,overshoot;
  UChar *x;
  clock_t end,start, end_real, start_real;
  struct tms r;
  double tot_time, tot_time_real;
  extern char *optarg;
  extern int optind, opterr, optopt;
  char *fnam, *sa_filename;
  char *lcp_filename,*bwt_filename;   // names for (optional) lcp and bwt files
  FILE *f;

  /* ------------ set default values ------------- */
  set_global_variables();
  print_sa=check_sa=0;
  sa_filename = NULL;
  lcp_filename = NULL;
  bwt_filename = NULL;

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  while ((c=getopt(argc, argv, "b:d:l:p:r:w:cvux:f:T:W:B:")) != -1) {
    switch (c) 
      {
      case 'b':
        bwt_filename = optarg; break;
      case 'c':
        check_sa++; break;
      case 'd':
        Anchor_dist = atoi(optarg); break;
      case 'l':
        Shallow_limit  = atoi(optarg); break;
      case 'p':
        print_sa = atoi(optarg); break;
      case 'x':
        _ds_Word_size = atoi(optarg); break;
      case 'v':
        _ds_Verbose++; break;
      case 'w':
        sa_filename = optarg; break;
      case 'f':
        Max_pseudo_anchor_offset = atoi(optarg); break;
      case 'r':
        B2g_ratio = atoi(optarg); break;
      case 'u':
        Update_anchor_ranks = 1; break;
      case 'T':
        Mk_qs_thresh = atoi(optarg); break;
      case 'W':
        lcp_filename = optarg; break;
      case 'B':
        Blind_sort_ratio = atoi(optarg); break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
      }
    num_opt++;
  }
  if(optind<argc)
    fnam=argv[optind];
  else {
    fprintf(stderr, "Usage:\n\t%s [-b bwtfile][-cuv][-d dist]",argv[0]);
    fprintf(stderr, "[-l len][-p num][-f maxoff][-r ratio]\n");
    fprintf(stderr, 
            "\t   [-T thresh][-w safile][-W lcpfile][-x wsize][-B ratio]");
    fprintf(stderr, " file\n\n");
    fprintf(stderr,"\t-b bwtfile  write bwt to bwtfile\n");    
    fprintf(stderr,
	    "\t-B ratio    blind_sort ratio [def. %d]\n",Blind_sort_ratio);
    fprintf(stderr,"\t-c          check the sa (could be very slow)\n");    
    fprintf(stderr,"\t-d dist     anchor distance [def. %d]\n",Anchor_dist);
    fprintf(stderr,"\t-f maxoff   Maximum offset for forward ");
    fprintf(stderr,"pseudo-anchors [def. %d]\n",Max_pseudo_anchor_offset);
    fprintf(stderr,
            "\t-l len      shallow sort limit [def. %d]\n",Shallow_limit);
    fprintf(stderr,
	    "\t-r ratio    bucket to group max ratio [def. %d]\n",B2g_ratio);
    fprintf(stderr,"\t-p num      print num char of each suffix [def. 0]\n");
    fprintf(stderr,
	    "\t-T thresh   Threshold for mk-qs [def. %d]\n", Mk_qs_thresh);
    fprintf(stderr,"\t-u          updates anchor ranks in get_rank()\n");
    fprintf(stderr,"\t-v          produces a verbose output\n");
    fprintf(stderr,"\t-w safile   write sa to safile\n");    
    fprintf(stderr,
            "\t-W lcpfile  check sa and write lcp to lcpfile (very slow)\n");
    fprintf(stderr,
	    "\t-x wsize    word size in mkqs (default %d)\n\n",_ds_Word_size); 
    return 0;
  }
  if(_ds_Verbose) {
    fprintf(stderr,"Command line: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n");
  }
  /* -------- check parameters ------------- */
  if(check_global_variables()) {
    exit(1);
  }

  /* ---------- open file and read text ----------- */
  if (! (f=fopen(fnam, "rb"))) {
    perror(fnam);
    return 1;
  }
  if (fseek(f, 0L, SEEK_END)) {
    perror(fnam);
    return 1;
  }
  n=ftell(f);
  if (n==0) {
    fprintf(stderr, "%s: file empty\n", fnam);
    return 0;
  }

  // ------ allocate memory for text and sa -------
  overshoot = compute_overshoot();
  p=malloc((n)*sizeof *p);               // sa
  x=malloc((n+overshoot)*sizeof *x);     // text
  if (! p || ! x) {
    fprintf(stderr, "malloc failed\n");
    return 1;
  }

  // ------------ read input text ---------------
  rewind(f); 
  c=fread(x, (size_t) 1, (size_t) n, f);
  // lseek(fileno(f),0,SEEK_SET); 
  // c=read(fileno(f), x, (size_t) n);
  if(c!=n) {
    fprintf(stderr,"Error in read() (%d vs %d) (main)\n",c,n);
    perror(fnam);
    return 1;
  }
  fclose(f);

  /* ---------  start measuring time ------------- */
  if(_ds_Verbose)
    fprintf(stderr,"Starting sa construction ... \n");
  start_real = times(&r);
  start  = (r.tms_utime+r.tms_stime);     /* user + system */
  ds_ssort(x, p, n);
  end_real = times(&r);
  end  = (r.tms_utime+r.tms_stime);     /* user + system */
 // tot_time =  ((double) (end-start))/CLK_TCK;
  //tot_time_real =  ((double) (end_real-start_real))/CLK_TCK;
  printf("Elapsed time: %.2f seconds (user+sys). Total real time: %.2f.\n", 
	 tot_time, tot_time_real);

  // --------------- write bwt to a file 
  if(bwt_filename!=NULL) 
    write_bwt(bwt_filename,x,p,n);

  // --------------- write sa to a file 
  if(sa_filename!=NULL) 
    write_sa(sa_filename,p,n);

  // --------------- write lcp to a file 
  if(lcp_filename!=NULL) 
    write_lcp(lcp_filename,x,p,n);

  // ------------ check sa --------
  if(check_sa) 
    check_sa_ordering(x,p,n,check_sa);

  // ----- display sa -------
  if(print_sa)   
    print_sa_onscreen(x,p,n,print_sa);

  // deallocate and exit
  free(x); free(p);
   return 0;
}


/* ********************************************************
   write bwt to filename: we output the last column of
   the bwt matrix (n bytes), followed by the position of input
   string (4 bytes)
   ******************************************************** */
void write_bwt(char *filename, UChar *text, int *sa, int text_size)
{
  FILE *bwt;
  Int32 i,eof_pos=0;

  if(_ds_Verbose)
    fprintf(stderr,"Writing bwt to file %s\n",filename);
  if((bwt=fopen(filename,"wb"))==NULL)
    perror(filename);

  putc(text[text_size-1],bwt);    // L[0] = Text[n-1]
  for(i=0;i<text_size;i++) 
    if(sa[i]!=0)
      putc(text[sa[i]-1],bwt);    // write char preceeding sa[i]
    else
      eof_pos = i;                // store position of EOF
  for(i=0;i<4;i++) {
    putc(eof_pos & 0xFF, bwt);    // write EOF
    eof_pos >>= 8;
  }
  fclose(bwt);
}

/* **********************************************************
   open filename and write p[0] .. p[n-1] using log(n) bits
   ********************************************************** */
void write_sa(char *filename, int *p, int n)
{
  int int_log2(int);
  void init_bit_buffer(void);
  void fbit_write(FILE *,int,int), fbit_flush( FILE * );
  FILE *sa;
  Int32 psize, i;

  if(_ds_Verbose)
    fprintf(stderr,"Writing sa to file %s\n",filename);
  if((sa=fopen(filename,"wb"))==NULL)
    perror(filename);

  init_bit_buffer();
  psize = int_log2(n);
  for(i=0;i<n;i++)
    fbit_write(sa,psize,p[i]);
  fbit_flush(sa);
  fclose(sa);
}


/* ********************************************************
   write lcp statistic to filename (plain ascii format)
   ******************************************************** */
void write_lcp(char *filename, UChar *x, int *p, int n)
{
  FILE *lcp;
  Int32 *stat, i, j, max_lcp=0, sum=0;
  unsigned long long sum_lcp=0;

  stat = (Int32 *) calloc(MAX_LCP_SIZE,sizeof(Int32)); // initialized to 0
  if(stat==NULL) {
    fprintf(stderr, "calloc failed (stat)\n");
    exit(1);
  }
  if(_ds_Verbose)
    fprintf(stderr,"Writing lcp stats to file %s\n",filename);
  if((lcp = fopen(filename,"w"))==NULL) 
     perror(filename);

  // computes lcp
  for(i=0;i<n-1;i++) {
    if (scmp3(x+p[i], x+p[i+1], & j, MIN(n-p[i], n-p[i+1]))>=0) {
      fprintf(stderr,"Error in sa file!\n");
      exit(1);
    }
    else {
      max_lcp = MAX(max_lcp,j);
      sum_lcp += j;
      if(j<MAX_LCP_SIZE)
	stat[j]++;   // one more lcp of length j
    }
  }
  // output lcp statistics
  fprintf(lcp,"Average lcp: %.2f\n",((double) sum_lcp)/(n-1));
  fprintf(lcp,"Maximum lcp: %d\n",max_lcp);
  if(max_lcp<MAX_LCP_SIZE) { 
    for(i=0;i<=max_lcp;i++) 
      if(stat[i]) {
	fprintf(lcp,"%10d %10d\n",i,stat[i]);
	sum += stat[i];
      }
    if(sum+1!=n) {
      fprintf(stderr,"Fatal error! Invalid lcp stats!\n");
      exit(1);
    }  
  }
  else {
    fprintf(stderr,"Unable to compute lcp stats. ");
    fprintf(stderr,"Please set MAX_LCP_SIZE to %d\n",max_lcp+1);
    exit(1);
  }
  fclose(lcp);
  free(stat);
}

// function for checking the sa (very slow) 
// if verbose>1 prints which suffixes are out of order
void check_sa_ordering(UChar *x, int *p, int n, int verbose)
{ 
  int i,j,wrong=0;

  printf("Checking...\n");
  for (i=0; i<n-1; ++i) {
    if (scmp3(x+p[i], x+p[i+1], & j, MIN(n-p[i], n-p[i+1]))>=0) {
      wrong++;
      if(verbose>1) {
	printf("---> i=%d  p[i]=%d  p[i+1]=%d\n", i, p[i], p[i+1]);
      }
    }  
  }
  if(wrong)
    printf("%d suffixes out of order!\n",wrong);
  else
    printf("done.\n");
}


// print the first print_sa char's of each suffix.
// if print_sa <0 print only the suffix starting position
void print_sa_onscreen(UChar *x, int *p, int n, int max_len)
{
  int i,j;

  for (i=0; i<n; ++i) {
    if(max_len<0)
      printf("%d\n", p[i]);
    else {
      printf("%3d] %3d \"", i, p[i]);
      for (j=p[i]; j<n && j-p[i]<max_len; ++j)
	pretty_putchar(x[j]);
      printf("\"\n");
    }
  }
}

// ----------------------------

/* *****************************************************
   Basic procedures to write bits using a Bit_buffer. 
   unwritten bits of Bit_buffer are the most significant ones.
   ***************************************************** */
typedef unsigned int uint32; 
uint32 Bit_buffer;       
int  Bit_buffer_size;  /* number of unread/unwritten bits in Bit_buffer */

// ***** Initializes the Bit_buffer to zero 
void init_bit_buffer(void)
{
  Bit_buffer= (uint32) 0;
  Bit_buffer_size=0;
}


// *****  Complete with zeroes the first byte of Bit_buffer 
//        This way, the content of Bit_buffer is entirely flushed out *****

void fbit_flush(FILE *f)
{
  void fbit_write24(FILE *f, int n, int vv);
  
  if(Bit_buffer_size!=0)
    fbit_write24(f, 8 - (Bit_buffer_size%8) , 0);  // pad with zero !
}


// ********* Write n (<= 24) bits taken from v. The content of 
//           Bit_buffer is flushed out until it contains <8 bits *********   
void fbit_write24(FILE *f, int n, int vv)  
{                              // v contains bits to read starting from
                               // the least significant bits
  uint32 v = (uint32) vv;

  assert(Bit_buffer_size<8);
  assert(n>0 && n<=24);
  assert( v < 1u <<(n+1) );
   
  /* ------- add n bits to Bit_buffer -------- */
  Bit_buffer_size += n;       // add first, to compute the correct shift
  Bit_buffer |= (v << (32 - Bit_buffer_size));  // compact to end of the buffer

  /* ------- flush Bit_buffer as much as possible ----- */
  while (Bit_buffer_size>=8) {            
    if( putc((Bit_buffer>>24),f) == EOF) {
      fprintf(stderr,"Error writing to output file -fbit_write-\n");
      exit(1);
    }
    Bit_buffer <<= 8;                       
    Bit_buffer_size -= 8;                 
  }                                           
} 

// ****** Write in file f the n bits taken from vv (possibly n > 24) 
void fbit_write(FILE *f, int n, int vv)
{  
  void fbit_write24(FILE *f,int n, int vv);  
  uint32 v = (uint32) vv;

  assert(n <= 32);
  if (n > 24){
    fbit_write24(f,n-24, (v>>24) & 0xffL);
    fbit_write24(f,24, v & 0xffffffL);
  } else {
    fbit_write24(f,n,v);
  }
}


int int_log2(int u)    // compute # bits to represent u
{
  int i = 1;
  int r = 1;
  
  while((i<=32) && (r<u)){
    r=2*r+1;
    i = i+1;
  }
    
  assert(i<=32);
  return i;
}






