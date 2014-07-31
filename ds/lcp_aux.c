/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   lcp_aux.c 
   Ver 1.0    28-apr-2004
   Functions for lightweight computation of the LCP array.

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
#include <assert.h>
#include "bwt_aux.h"

/* ==================================================================
   Description of the data:
   text: [0,n-1] -> [0,255]
   suffix array: [1,n] -> [0,n-1], 
      sa[i] is the starting position of the i-th suffix 
      in lexicographic order
   lcp: [1,n]->[0,n-1],
      lcp[i] is the length of the LCP between t[sa[i],n-1] and 
      its predecessor in the lexicographic order (which is t[sa[i-1],n-1]).
      lcp[1] is undefined since t[sa[1],n-1] has no predecessor
  ======================================================================== */

/* ***********************************************************************
    Compute the lcp array from the suffix array
    using the algorithm by Kasai et al. (CPM '01)
    input
      t[0,n-1] input text 
      sa[1,n] suffix array
    return
      lcp array, or NULL if an error (e.g. out of memory) occurs.
   *********************************************************************** */
int *_lcp_sa2lcp_13n(uint8 *t, int n, int *sa)
{
  int i,h,k,j, *lcp, *rank;

  lcp = (int*)malloc((n+1)*sizeof(int));
  rank = (int*)malloc(n*sizeof(int));
  if(lcp==NULL || rank==NULL)
    return NULL;

  // compute rank = sa^{-1}
  for(i=1;i<=n;i++) rank[sa[i]] = i;
  // traverse suffixes in rank order
  h=0;
  for(i=0;i<n;i++) {
    k = rank[i];          // rank of s[i ... n-1]
    if(k>1) {
      j = sa[k-1];        // predecessor of s[i ... n-1]
      while(i+h<n && j+h<n && t[i+h]==t[j+h])
	h++;
      lcp[k] = h;
    }
    if(h>0) h--;
  }
  free(rank);
  return lcp;
}

/* ***********************************************************************
   space economical computation of the lcp array
   This procedure is similar to the previous one, but we compute the rank
   of i using the rank_next map (stored in the rn array) instead of
   the rank array. The advantage is that as soon as we have read rn[k] that
   position is no longer needed and we can use it to store the lcp.
   Thus, rn[] and lcp[] share the same memory and the overall space
   requirement of the procedure is 9n instead of 13n.
    input
      t[0,n-1] input text 
      sa[1,n] suffix array
      occ[0,_BW_ALPHA_SIZE-1]  # of occurrences of each char
    return
      lcp array, or NULL if an error (e.g. out of memory) occurs.
    space
     9n
   *********************************************************************** */
int *_lcp_sa2lcp_9n(uint8 *t, int n, int *sa, int *occ)
{
  int i,h,j,k,nextk=-1;
  int *rn, *lcp;

  lcp = (int*)malloc((n+1)*sizeof(int));
  if(lcp==NULL)
    return NULL;  // out of memory

  // rn and lcp are pointers to the same array
  rn = lcp;
  // compute rank_next map
  k = _bw_sa2ranknext(t, n, sa, occ, rn); // k is the rank of s
  for(h=i=0; i<n; i++,k=nextk) {
    assert(k>0 && i==sa[k]); 
    nextk=rn[k];          // read nextk before it is overwritten
    if(k>1) {
      // compute lcp between suffixes of rank k and k-1 (recall i==sa[k])
      j = sa[k-1];
      while(i+h<n && j+h<n && t[i+h]==t[j+h])
	h++;
      lcp[k]=h;           // h is the lcp, write it to lcp[k];  
    }
    if(h>0) h--;
  }
  assert(nextk==0);        // we have reached the last suffix s[n-1]
  return lcp;
}



/* ***********************************************************************
   space economical computation of the lcp array based on the bwt
    input
      t[0,n-1] input text 
      b->bwt[0,n] b->eof_pos   bwt of t[0,n-1]
      sa[1,n] suffix array
      occ[0,_BW_ALPHA_SIZE-1]  # of occurrences of each char
    output 
      the lcp array is overwritten to sa
    return
      num = number of extra int32 used by the algorithm, 
      or -1 if an error (e.g. out of memory) occurs.
    space
     (6 + delta)n, where delta=(4*num)/n
   *********************************************************************** */
int _lcp_sa2lcp_6n(uint8 *t, bwt_data *b, int *sa, int *occ)
{
  int i,h,j,k,num,num2,nextk=-1,n=b->size;
  int *rn=sa, *lcp=sa, *sa_aux;

  // ----- dirty trick: make b->bwt[eof_pos] different from 
  //       b->bwt[eof_pos-1] and b->bwt[eof_pos+1]
  assert(b->eof_pos!=0);
  b->bwt[b->eof_pos] = (uint8) (b->bwt[b->eof_pos-1] ^ 1); // change 0th bit
  if(b->eof_pos<n && b->bwt[b->eof_pos]==b->bwt[b->eof_pos+1])
    b->bwt[b->eof_pos] = (uint8) (b->bwt[b->eof_pos] ^ 2); // change 1st bit
    
  // ----- count # of sa positions that we need later
  for(num=0,i=2; i<=n; i++) 
    if(b->bwt[i-1]!=b->bwt[i]) num++;
  assert(num>0 && num<=n);

  // ---------- alloc extra memory -----------
  sa_aux = (int *) malloc(num*sizeof(int));
  if(sa_aux==NULL)
    return -1;
 
  // ---------- compute the ranknext map from the bwt ----
  k = _bw_bwt2ranknext(b, occ, rn);   // k is the rank of t[0,n-1]

  // ---------- convert rn->sa (in place) and save useful sa positions
  for(num2=i=0; k!=0; ) {
    assert(i<n);
    if( k>1 && (b->bwt[k-1]!=b->bwt[k]) )
      sa_aux[num2++] = k-1;
    nextk = rn[k];
    sa[k] = i++;
    k = nextk;
  }
  assert(num2==num);      // all needed positions have been recorded
  assert(i==n);           // we have reached the last suffix t[n-1]

  // ------ save sa values in sa_aux
  for(i=0;i<num;i++) {
    assert(sa_aux[i]>0 && sa_aux[i]<=n);
    sa_aux[i]=sa[sa_aux[i]];
  }

  // ------ now compute the lcp
  k = _bw_bwt2ranknext(b, occ, rn);   // k is the rank of t[0,n-1]
  for(num2=h=i=0; i<n; i++,k=nextk) {
    assert(k>0);
    nextk=rn[k];               // read nextk before it is overwritten
    if(k>1) {
      if(b->bwt[k]!=b->bwt[k-1]) {
	j = sa_aux[num2++];    // retrieve sa[k-1]
	while(i+h<n && j+h<n && t[i+h]==t[j+h]) 
          h++;                 // extend lcp if possible
      }
      lcp[k]=h;                // h is the lcp, write it to lcp[k];  
    }
    if(h>0) h--;
  }
  assert(nextk==0);        // we have reached the last suffix t[n-1]
  assert(num2==num);
  free(sa_aux);
  return num;
}


/* ***********************************************************************
    Compute the lcp array from the suffix array
    using the algorithm by V. Makkinen ( '01)
    input
      t[0,n-1] input text 
      sa[1,n] suffix array
    return
      lcp array, or NULL if an error (e.g. out of memory) occurs.
   *********************************************************************** */
#define MARKER (1<<31)
int *_lcp_sa2lcp_9125n(uint8 *t, int n, int *sa)
{
  int i,h,k,j, *lcp, *rank;
  int nextj,lcpi,nextlcp; 

  lcp = (int*)malloc((n+1)*sizeof(int));
  if(lcp==NULL)
    return NULL;

  // --- in the first part of the procedure rank[] and lcp[] coincide
  rank=lcp;
  // compute rank = sa^{-1}
  for(i=1;i<=n;i++) rank[sa[i]] = i;
  // traverse suffixes in rank order
  for(h=i=0;i<n;i++) {
    k = rank[i];          // rank of s[i ... n-1]
    if(k>1) {
      j = sa[k-1];        // predecessor of s[i ... n-1]
      while(i+h<n && j+h<n && t[i+h]==t[j+h])
	h++;
      lcp[i] = h;         // h is lcp[k]=lcp[rank[i]]
    }
    if(h>0) h--;
  }
  // now the lcp values have been computed but they are in
  // in the wrong positions....
  rank = sa;  // from now on rank[] and lcp[] coincide

  /* ---------------------------------------------------------------------
     Now we transform (in place) sa[] values into the rank values.
     We use to MSB to mark the entries which have been already overwritten
     --------------------------------------------------------------------- */
  sa[0]=0;                               // mark this entry as not written
  for(k=n;k>=1;k--) {
    if(!(sa[k] & MARKER)) {              // if sa[k] not marked enter loop
      i=k; j=sa[i];              
      while(i>0 && !(sa[j] & MARKER)) {  // loop invariant: j = sa[i]
        assert(j<n);
        nextj=sa[j];                     // save content of sa[j]
        sa[j]= (i | MARKER);             // rank[j] gets i and is marked
        i=j; j=nextj;                    // j=nextj'=sa[j']=sa[i] 
      }
    }
  }

  /* ------------------------------------------------------------------- 
     now we transform (in place) rank values to sa values and 
     move lcp values in their proper location.
     Note that the entries to be overwritten are now those with the MSB set
     ------------------------------------------------------------------- */
  rank[n]=MARKER;                // mark this entry as not written
  for(k=0;k<n;k++) {
    if((rank[k] & MARKER)) {     // if rank[k] is marked enter loop
      i=k; j=rank[i] & ~MARKER;
      lcpi = lcp[i];
      while(i<n && (rank[j] & MARKER)) {
        // loop invariant: j = rk[i], lcpi=lcp[i]=truelcp(j)
        assert(j>0);
        nextj=rank[j] & ~MARKER; // save content of rank[j]
        rank[j]=i;               // rank[j] gets i and is unmarked
        nextlcp=lcp[j];          // save lcp[j]
        lcp[j]=lcpi;             // write truelcp(j) in the right place
        i=j; j=nextj;            // j=nextj'=rank[j']=rank[i]
        lcpi = nextlcp;          // lcpi=lcp[j']=lcp[i]
      }
    }
  }
  // ---- now lcp and sa values are in the proper order
  return lcp;
}



