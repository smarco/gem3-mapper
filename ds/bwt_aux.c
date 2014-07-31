/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt_aux.c 
   Ver 1.0    3-mar-2004
   Functions for computing and transforming various representations 
   of text and suffix array: t/sa/rank_next/rank_prev/bwt

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
#include <assert.h>
#include "bwt_aux.h"

/* ==================================================================
   Description of the data:
   text: [0,n-1] -> [0,255]
   suffix array: [1,n] -> [0,n-1], 
      sa[i] is the starting position of the i-th suffix 
      in lexicographic order
   rank: [0,n-1] -> [1,n]
      inverse of sa.  rank[i] is rank of t[i,n-1] in the lexicographic order.
      rank is never used explicitly in the following routines, but
      is is useful for the definition of rank_next and rank_prev
   rank_next: [1,n] -> [0,n]
      rnext[i] = rank[sa[i]+1]; if sa[i]=n-1 (i.e. i=rank[n-1]) rnext[i]=0 
      note that rnext[i]!=rank[0] for every i.
   rank_prev: [1,n] -> [0,n]
      rprev[i] = rank[sa[i]-1]; if sa[i]=0 (i.e. i=rank[0]) -> rprev[i]=0 
      note that rprev[i]!=rank[n-1] for every i.
   The BWT is represented by a three element struct:
     bwt: [0,n] -> [0,255], 
          bwt[0]=t[n-1], bwt[rank[0]] undefined, bwt[i]=t[sa[i]-1] 
     size = n;
     eof_pos = rank[0];
  ======================================================================== */
 
    
/* *******************************************************************
   "Standard" computation of the BWT given the suffix array.
   input
     n size of text
     t[0,n-1] input text
     sa[1,n] suffix array
     bwt_data *b: at call b->bwt should have been already allocated
   output
     bwt_data *b
   space: 6n
   ******************************************************************* */   
void _bw_sa2bwt(uchar *t, int32 n, int32 *sa, bwt_data *b)
{
  int i;

  assert(b->bwt!=NULL);
  b->size = n;
  b->bwt[0] = t[n-1];
  for(i=1;i<=n;i++)
    if(sa[i]>0)
      b->bwt[i]=t[sa[i]-1];
    else
      b->eof_pos=i;
}


/* =======================================================================
   functions based on the rank_next
   ======================================================================= */

/* *******************************************************************
   input 
     bwt_data, occ
   output
     rank_next
   return 
     r0 = rank of t[0,n-1], sa[r0]=0 
   space 5n
   ******************************************************************* */
int32 _bw_bwt2ranknext(bwt_data *b, int32* occ, int32 *rank_next)
{
  int32 i,j,n,c,r0=0;
  int32 count[_BW_ALPHA_SIZE];

  // --- occ -> count
  count[0]=0;
  for(i=1;i<_BW_ALPHA_SIZE;i++)
    count[i] = count[i-1] + occ[i-1];
  // --- bwt -> rank_next
  n = b->size;
  for(i=0;i<=n;i++) 
    if(i == b->eof_pos)
      r0 = i;
    else {
      c = b->bwt[i];
      j = ++count[c];
      rank_next[j]=i;
    }
  assert(r0>0);
  return r0;
}

/* *******************************************************************
  input 
     n size of text
     t[0,n-1] input text
     sa[1,n] suffix array
     occ
  output
    rank_next
  return 
    r0 = rank of t[0,n-1], sa[r0]=0 
  space 8n
   ******************************************************************* */
int32 _bw_sa2ranknext(uchar *t, int32 n, int32 *sa, int32 *occ, 
                      int32 *rank_next)
{
  int32 i,j,c,r0=0,count[_BW_ALPHA_SIZE];

  // --- occ -> count
  count[0]=0;
  for(i=1;i<_BW_ALPHA_SIZE;i++)
    count[i] = count[i-1] + occ[i-1];
  // --- sa+t -> rank_next
  j = ++count[t[n-1]];       // this is bwt[0]
  rank_next[j]=0;
  for(i=1;i<=n;i++) 
    if(sa[i] == 0)
      r0 = i;
    else {
      c = t[sa[i]-1];
      j = ++count[c];
      rank_next[j]=i;
    }
  assert(r0>0);
  return r0;
}


/* *******************************************************************
   input 
     rank_next, r0, bwt_data
   output
     t[]
   space 6n
   ******************************************************************* */
void _bw_ranknext2t(int32 *rank_next, int32 r0, bwt_data *b, uchar *t)
{ 
  int32 k,i;

  k = r0; i=0;
  do {
    k = rank_next[k];
    t[i++] = b->bwt[k];
  } while(k!=0);
  assert(i==b->size);
}


/* *******************************************************************
   input 
     rank_next, r0
   output
     sa
   space 4n (in place)/8n
   ******************************************************************* */
void _bw_ranknext2sa(int32 *rank_next, int32 r0, int32 *sa)
{
  int32 k, nextk, i;

  // --- if sa==NULL sa entries are written to rank_next 
  if(sa==NULL)
    sa=rank_next;
  // --- compute sa
  k=r0; i=0;
  while(k!=0) {
    nextk = rank_next[k];
    sa[k] = i++;
    k = nextk;
  }
  // assert(i==n); we cannot test this since n is not passed to the procedure
}




/* ===================================================================
   functions base on rank_prev
   =================================================================== */

/* *******************************************************************
   input 
     bwt_data, occ
   output
     rank_prev
   return 
     rn1 = rank of t[n-1], sa[rn1]=n-1 
   space 5n
   ******************************************************************* */
int32 _bw_bwt2rankprev(bwt_data *b, int32* occ, int32 *rank_prev)
{
  int32 i,j,n,c,rn1;
  int32 count[_BW_ALPHA_SIZE];

  // --- occ -> count
  count[0]=0;
  for(i=1;i<_BW_ALPHA_SIZE;i++)
    count[i] = count[i-1] + occ[i-1];
  // --- bwt -> rank_next
  n = b->size;
  rn1 = ++count[b->bwt[0]];
  for(i=1;i<=n;i++) 
    if(i == b->eof_pos)
      rank_prev[i]=0;
    else {
      c = b->bwt[i];
      j = ++count[c];
      rank_prev[i]=j;
    }
  return rn1;
}

/* *******************************************************************
   input 
     t, sa
   output
     rank_prev
   return 
     rn1 = rank of t[n-1], sa[rn1]=n-1 
   space 8n/ 4n (in place)
   ******************************************************************* */
int32 _bw_sa2rankprev(uchar *t, int32 n, int32 *sa, int32 *occ, 
                      int32 *rank_prev)
{
  int32 i,j,c,rn1;
  int32 count[_BW_ALPHA_SIZE];

  // --- if rank_prev==NULL rank_prev entries are overwritten in sa
  if(rank_prev==NULL)
    rank_prev=sa;
  // --- occ -> count
  count[0]=0;
  for(i=1;i<_BW_ALPHA_SIZE;i++)
    count[i] = count[i-1] + occ[i-1];
  // --- (t+sa) -> rank_prev
  j = ++count[t[n-1]];
  rn1 = j;
  for(i=1;i<=n;i++) 
    if(sa[i] == 0)
      rank_prev[i]=0;
    else {
      c = t[sa[i]-1];     // this is bwt[i];
      j = ++count[c];
      rank_prev[i]=j;
    }
  return rn1;
}

/* *******************************************************************
   input 
     rank_prev, rn1, bwt_data
   output
     t[]
   space 6n
   ******************************************************************* */
void _bw_rankprev2t(int32 *rank_prev, int32 rn1, bwt_data *b, uchar *t)
{ 
  int32 k,i,n;

  n = b->size;          // size of t     
  t[n-1] = b->bwt[0];   // last char of t is bwt[0] 
  k = rn1; i=n-1;       // init loop
  do {
    assert(k!=0);
    t[--i] = b->bwt[k]; // reconstruct t backward
    k = rank_prev[k];
  } while(i>0);
  assert(k==b->eof_pos);
  assert(rank_prev[k]==0);
}

/* *******************************************************************
   input 
     rank_prev, n, rn1 (rank[n-1])
   output
     sa
   space 4n (in place)/8n
   ******************************************************************* */
void _bw_rankprev2sa(int32 *rank_prev, int32 n, int32 rn1, int32 *sa)
{
  int32 k, prevk, i;

  // --- if sa==NULL sa entries are written to rank_prev 
  if(sa==NULL)
    sa=rank_prev;
  // --- compute sa
  k=rn1; i=n-1;
  while(k!=0) {
    prevk = rank_prev[k];
    sa[k] = i--;
    k = prevk;
  }
  assert(i== -1);
}


/* *******************************************************************
   input 
     rank_prev, rn1 (rank[n-1])
   output
     rank_next
   return 
     r0 (rank[0])
   space: 4n (in place)/8n
   ******************************************************************* */
int _bw_rprev2rnext(int32 *rank_prev, int32 rn1, int32 *rank_next)
{
  int32 rprec, rcur,rnext;

  // --- if rank_next==NULL rank_next entries are written to rank_prev 
  if(rank_next==NULL)
    rank_next=rank_prev;

  rnext = rn1;
  rcur = rank_prev[rnext];
  // inside the loop we have the invariant rnext=rank_next(rcur)  
  while(rcur!=0) {
    rprec = rank_prev[rcur];    // save rank_prev[rcur]
    rank_next[rcur] = rnext;    // true by the invariant
    rnext = rcur;               // update 
    rcur = rprec;               // the invariant still hold
  }
  rank_next[rn1]=0;             
  return rnext;                 // since rcur=rank_prev[rnext]
                                // rcur==0 implies rnext=rank[0] 
}

