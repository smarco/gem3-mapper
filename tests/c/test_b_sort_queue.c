/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2026 by Jacob Schreiber  <jacob.schreiber@umassmed.edu>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Jacob Schreiber <jacob.schreiber@umassmed.edu> (unit tests, written with Claude Code)
 * DESCRIPTION:
 *   The sort template (include/utils/vector_sort.h: selection sort,
 *   heapsort, quicksort with a bounded call stack), the sorts gem-mapper
 *   instantiates from it for candidates and alignment regions, the priority
 *   queue, and the alignment-region accessors, against sorted copies made
 *   with qsort and hand-worked values. Sorts are checked for order and for
 *   being a permutation of the input; the relative order of equal keys is
 *   not specified and not checked.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     vector_sort_filtering_positions(vector_t*)  (filtering_position_t*, by text_end_position)
 *     vector_sort_align_distance(vector_t*)       (filtering_region_t*, by distance_min_bound, then text_begin_position)
 *     vector_sort_scaffold_coverage(vector_t*)    (filtering_region_t*, by scaffolding_coverage, descending)
 *     buffer_sort_match_alignment_region(match_alignment_region_t*,int)  (by text_begin)
 *     pqueue_new(uint64_t), pqueue_push_(pqueue_t*,void*,uint64_t), pqueue_pop_(pqueue_t*),
 *     pqueue_top_priority(pqueue_t*), pqueue_delete(pqueue_t*)
 *     match_alignment_region_{init,get_*,set_*,text_coverage,text_distance,text_overlap,
 *         cmp_key_position,cmp_text_position}
 */

#include "gem3_test.h"
#include "utils/priority_queue.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/region/filtering_region.h"
#include "matches/align/match_alignment_region.h"

// Instances defined in filtering_candidates.c and match_scaffold.c (no header declares them)
void vector_sort_filtering_positions(vector_t* const vector);
void vector_sort_align_distance(vector_t* const vector);
void vector_sort_scaffold_coverage(vector_t* const vector);
void buffer_sort_match_alignment_region(match_alignment_region_t* const vector_mem,const int num_positions);

// The template on plain integers, instantiated here
#define VECTOR_SORT_NAME                 gt_int64
#define VECTOR_SORT_TYPE                 int64_t
#define VECTOR_SORT_CMP(a,b)             ((*(a) > *(b)) - (*(a) < *(b)))
#include "utils/vector_sort.h"

#define MAX_N 5000

static int cmp_int64(const void* a,const void* b) {
  const int64_t x = *(const int64_t*)a, y = *(const int64_t*)b;
  return (x > y) - (x < y);
}

// Input shapes: random, few distinct values, sorted, reversed, constant, organ pipe
static void fill(gt_rng_t* const rng,int64_t* const v,const int n,const int shape) {
  int i;
  for (i=0;i<n;++i) {
    switch (shape) {
      case 0: v[i] = (int64_t)gt_rng_below(rng,1000000) - 500000; break;
      case 1: v[i] = (int64_t)gt_rng_below(rng,4); break;
      case 2: v[i] = i; break;
      case 3: v[i] = n-i; break;
      case 4: v[i] = 7; break;
      default: v[i] = (i < n/2) ? i : n-i; break;
    }
  }
}

static const int sizes[] = { 0, 1, 2, 3, 10, 11, 12, 13, 31, 100, 257, 1000, 4999 };
#define NUM_SIZES (int)(sizeof(sizes)/sizeof(sizes[0]))

/*
 * The template: buffer_sort, vector_sort and heapsort agree with qsort on
 * every shape and size, including the sorted and constant inputs that drive
 * the quicksort call stack to its depth limit (heapsort fallback).
 */
GT_TEST(sort_template_matches_qsort) {
  gt_rng_t rng = gt_rng_new(11);
  static int64_t v[MAX_N], ref[MAX_N];
  vector_t* const vec = vector_new(MAX_N,int64_t);
  int s, shape;
  for (s=0;s<NUM_SIZES;++s) {
    const int n = sizes[s];
    for (shape=0;shape<6;++shape) {
      fill(&rng,v,n,shape);
      memcpy(ref,v,n*sizeof(int64_t));
      qsort(ref,n,sizeof(int64_t),cmp_int64);
      // buffer_sort
      int64_t* const a = vector_get_mem(vec,int64_t);
      memcpy(a,v,n*sizeof(int64_t));
      buffer_sort_gt_int64(a,n);
      GT_CHECK_MSG(n==0 || memcmp(a,ref,n*sizeof(int64_t))==0,"buffer_sort n=%d shape=%d",n,shape);
      // vector_sort
      memcpy(a,v,n*sizeof(int64_t));
      vector_set_used(vec,n);
      vector_sort_gt_int64(vec);
      GT_CHECK_MSG(n==0 || memcmp(vector_get_mem(vec,int64_t),ref,n*sizeof(int64_t))==0,"vector_sort n=%d shape=%d",n,shape);
      GT_CHECK_EQ_UINT(vector_get_used(vec),n);
      // heapsort and selection alone
      if (n > 0) {
        memcpy(a,v,n*sizeof(int64_t));
        vector_sort_heapsort_gt_int64(a,n);
        GT_CHECK_MSG(memcmp(a,ref,n*sizeof(int64_t))==0,"heapsort n=%d shape=%d",n,shape);
      }
      if (n <= 257) {
        memcpy(a,v,n*sizeof(int64_t));
        vector_sort_selection_gt_int64(a,n);
        GT_CHECK_MSG(n==0 || memcmp(a,ref,n*sizeof(int64_t))==0,"selection n=%d shape=%d",n,shape);
      }
    }
  }
  vector_delete(vec);
}

/*
 * gem-mapper's instances sort by their documented keys and permute their
 * input (each element present once).
 */
GT_TEST(sort_instances_order_and_permutation) {
  gt_rng_t rng = gt_rng_new(22);
  static filtering_position_t positions[MAX_N];
  static filtering_region_t regions[1200];
  static match_alignment_region_t aregions[MAX_N];
  static uint8_t seen[MAX_N];
  vector_t* const vec = vector_new(MAX_N,void*);
  int s, shape, i;
  for (s=0;s<NUM_SIZES;++s) {
    const int n = sizes[s];
    for (shape=0;shape<6;++shape) {
      static int64_t keys[MAX_N], keys2[MAX_N];
      fill(&rng,keys,n,shape);
      fill(&rng,keys2,n,(shape+1)%6);
      // filtering positions by text_end_position
      filtering_position_t** const pp = vector_get_mem(vec,filtering_position_t*);
      for (i=0;i<n;++i) {
        positions[i].text_end_position = (uint64_t)(keys[i]+1000000);
        pp[i] = positions+i;
      }
      vector_set_used(vec,n);
      vector_sort_filtering_positions(vec);
      memset(seen,0,n);
      bool ok = true;
      for (i=0;i<n;++i) {
        const int idx = (int)(pp[i]-positions);
        ok = ok && idx>=0 && idx<n && !seen[idx];
        if (idx>=0 && idx<n) seen[idx] = 1;
        if (i>0) ok = ok && pp[i-1]->text_end_position <= pp[i]->text_end_position;
      }
      GT_CHECK_MSG(ok,"filtering_positions n=%d shape=%d",n,shape);
      // filtering regions (at most 1200 of them)
      const int nr = MIN(n,1200);
      filtering_region_t** const rp = vector_get_mem(vec,filtering_region_t*);
      for (i=0;i<nr;++i) {
        regions[i].alignment.distance_min_bound = (uint64_t)(keys[i] % 7 + 7) % 7;
        regions[i].text_begin_position = (uint64_t)(keys2[i]+1000000);
        regions[i].match_scaffold.scaffolding_coverage = (uint64_t)(keys[i]+1000000);
        rp[i] = regions+i;
      }
      vector_set_used(vec,nr);
      vector_sort_align_distance(vec);
      memset(seen,0,nr);
      ok = true;
      for (i=0;i<nr;++i) {
        const int idx = (int)(rp[i]-regions);
        ok = ok && idx>=0 && idx<nr && !seen[idx];
        if (idx>=0 && idx<nr) seen[idx] = 1;
        if (i>0) {
          const filtering_region_t* const a = rp[i-1];
          const filtering_region_t* const b = rp[i];
          ok = ok && (a->alignment.distance_min_bound < b->alignment.distance_min_bound ||
              (a->alignment.distance_min_bound == b->alignment.distance_min_bound &&
               a->text_begin_position <= b->text_begin_position));
        }
      }
      GT_CHECK_MSG(ok,"align_distance n=%d shape=%d",nr,shape);
      for (i=0;i<nr;++i) rp[i] = regions+i;
      vector_set_used(vec,nr);
      vector_sort_scaffold_coverage(vec);
      memset(seen,0,nr);
      ok = true;
      for (i=0;i<nr;++i) {
        const int idx = (int)(rp[i]-regions);
        ok = ok && idx>=0 && idx<nr && !seen[idx];
        if (idx>=0 && idx<nr) seen[idx] = 1;
        if (i>0) ok = ok && rp[i-1]->match_scaffold.scaffolding_coverage >= rp[i]->match_scaffold.scaffolding_coverage;
      }
      GT_CHECK_MSG(ok,"scaffold_coverage n=%d shape=%d",nr,shape);
      // alignment regions by text_begin (sorted in place)
      for (i=0;i<n;++i) {
        match_alignment_region_init(aregions+i,match_alignment_region_exact,0,i,1,
            i,i+1,(uint64_t)(keys[i]+1000000),(uint64_t)(keys[i]+1000010));
      }
      buffer_sort_match_alignment_region(aregions,n);
      memset(seen,0,n);
      ok = true;
      for (i=0;i<n;++i) {
        const uint64_t id = match_alignment_region_get_cigar_buffer_offset(aregions+i);
        ok = ok && id<(uint64_t)n && !seen[id];
        if (id<(uint64_t)n) seen[id] = 1;
        ok = ok && match_alignment_region_get_text_begin(aregions+i)==(uint64_t)(keys[id]+1000000);
        if (i>0) ok = ok && match_alignment_region_cmp_text_position(aregions+i-1,aregions+i) <= 0;
      }
      GT_CHECK_MSG(ok,"match_alignment_region n=%d shape=%d",n,shape);
    }
  }
  vector_delete(vec);
}

/*
 * Priority queue: pops come out in ascending priority; the top is the
 * minimum; an empty queue pops NULL with top priority UINT64_MAX.
 */
GT_TEST(pqueue_orders_by_priority) {
  gt_rng_t rng = gt_rng_new(33);
  static int64_t items[MAX_N];
  static uint64_t prio[MAX_N], sorted[MAX_N];
  int round;
  for (round=0;round<40;++round) {
    pqueue_t* const pq = pqueue_new(round % 3 == 0 ? 0 : 64);
    GT_CHECK(pqueue_is_empty(pq));
    GT_CHECK_EQ_UINT(pqueue_top_priority(pq),UINT64_MAX);
    GT_CHECK(pqueue_pop(pq,int64_t)==NULL);
    const int n = (int)gt_rng_range(&rng,1,round < 20 ? 30 : 3000);
    int i;
    for (i=0;i<n;++i) {
      prio[i] = (round % 4 == 0) ? gt_rng_below(&rng,5) : gt_rng_below(&rng,100000);
      items[i] = i;
      sorted[i] = prio[i];
      int64_t* const item = items+i; // pqueue_push casts its argument unparenthesized
      pqueue_push(pq,item,prio[i]);
    }
    GT_CHECK_EQ_UINT(pqueue_get_num_elements(pq),n);
    qsort(sorted,n,sizeof(uint64_t),cmp_int64);
    GT_CHECK_EQ_UINT(pqueue_top_priority(pq),sorted[0]);
    bool ok = true;
    for (i=0;i<n;++i) {
      const uint64_t top = pqueue_top_priority(pq);
      int64_t* const e = pqueue_pop(pq,int64_t);
      ok = ok && e!=NULL && top==sorted[i] && prio[*e]==sorted[i];
    }
    GT_CHECK_MSG(ok,"round %d: pops out of priority order",round);
    GT_CHECK(pqueue_is_empty(pq));
    GT_CHECK(pqueue_pop(pq,int64_t)==NULL);
    pqueue_delete(pq);
  }
  // Interleaved pushes and pops against a sorted reference
  pqueue_t* const pq = pqueue_new(0);
  static uint64_t ref[MAX_N];
  int nref = 0, op;
  bool ok = true;
  for (op=0;op<4000;++op) {
    if (nref==0 || gt_rng_chance(&rng,3,5)) {
      const uint64_t p = gt_rng_below(&rng,1000);
      int64_t* const item = items+(op % MAX_N);
      *item = (int64_t)p;
      pqueue_push(pq,item,p);
      ref[nref++] = p;
    } else {
      qsort(ref,nref,sizeof(uint64_t),cmp_int64);
      int64_t* const e = pqueue_pop(pq,int64_t);
      ok = ok && e!=NULL && (uint64_t)*e==ref[0];
      memmove(ref,ref+1,(nref-1)*sizeof(uint64_t));
      --nref;
    }
  }
  GT_CHECK_MSG(ok,"interleaved push/pop out of order");
  GT_CHECK_EQ_UINT(pqueue_get_num_elements(pq),nref);
  pqueue_delete(pq);
}

/*
 * Alignment-region accessors and text arithmetic
 */
GT_TEST(match_alignment_region_accessors) {
  match_alignment_region_t a, b;
  match_alignment_region_init(&a,match_alignment_region_approximate,3,40,5,10,60,1000,1052);
  GT_CHECK_EQ_INT(match_alignment_region_get_type(&a),match_alignment_region_approximate);
  GT_CHECK_EQ_UINT(match_alignment_region_get_error(&a),3);
  GT_CHECK_EQ_UINT(match_alignment_region_get_cigar_buffer_offset(&a),40);
  GT_CHECK_EQ_UINT(match_alignment_region_get_cigar_length(&a),5);
  GT_CHECK_EQ_UINT(match_alignment_region_get_key_begin(&a),10);
  GT_CHECK_EQ_UINT(match_alignment_region_get_key_end(&a),60);
  GT_CHECK_EQ_UINT(match_alignment_region_get_text_begin(&a),1000);
  GT_CHECK_EQ_UINT(match_alignment_region_get_text_end(&a),1052);
  GT_CHECK_EQ_INT(match_alignment_region_text_coverage(&a),52);
  match_alignment_region_set_type(&a,match_alignment_region_exact);
  GT_CHECK_EQ_INT(match_alignment_region_get_type(&a),match_alignment_region_exact);
  match_alignment_region_set_cigar_buffer_offset(&a,7);
  match_alignment_region_set_cigar_length(&a,2);
  match_alignment_region_set_key_begin(&a,11);
  match_alignment_region_set_key_end(&a,59);
  match_alignment_region_set_text_begin(&a,1001);
  match_alignment_region_set_text_end(&a,1049);
  GT_CHECK_EQ_UINT(match_alignment_region_get_cigar_buffer_offset(&a),7);
  GT_CHECK_EQ_UINT(match_alignment_region_get_cigar_length(&a),2);
  GT_CHECK_EQ_UINT(match_alignment_region_get_key_begin(&a),11);
  GT_CHECK_EQ_UINT(match_alignment_region_get_key_end(&a),59);
  GT_CHECK_EQ_UINT(match_alignment_region_get_text_begin(&a),1001);
  GT_CHECK_EQ_UINT(match_alignment_region_get_text_end(&a),1049);
  // b after a with a gap of 21 text bases; then overlapping; then touching
  match_alignment_region_init(&b,match_alignment_region_exact,0,0,1,70,90,1070,1090);
  GT_CHECK_EQ_INT(match_alignment_region_text_distance(&a,&b),21);
  GT_CHECK(!match_alignment_region_text_overlap(&a,&b));
  GT_CHECK(match_alignment_region_cmp_text_position(&a,&b) < 0);
  GT_CHECK(match_alignment_region_cmp_key_position(&b,&a) > 0);
  match_alignment_region_set_text_begin(&b,1048);
  GT_CHECK(match_alignment_region_text_overlap(&a,&b));
  GT_CHECK(match_alignment_region_text_overlap(&b,&a));
  GT_CHECK_EQ_INT(match_alignment_region_text_distance(&a,&b),-1);
  match_alignment_region_set_text_begin(&b,1049);
  GT_CHECK(!match_alignment_region_text_overlap(&a,&b));
  GT_CHECK_EQ_INT(match_alignment_region_text_distance(&a,&b),0);
  GT_CHECK_EQ_INT(match_alignment_region_cmp_text_position(&a,&a),0);
}

GT_CASES(
  GT_CASE(sort_template_matches_qsort),
  GT_CASE(sort_instances_order_and_permutation),
  GT_CASE(pqueue_orders_by_priority),
  GT_CASE(match_alignment_region_accessors)
)
