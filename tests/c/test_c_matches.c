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
 *   Single-end match handling: CIGAR vectors and their distances, match
 *   counters (strata), match metrics, the match container and its sort
 *   orders, match classification, MAPQ predictors, match selection
 *   (-m/-M) and the single-end MAPQ models.
 *
 *   GEM CIGAR elements: cigar_match (M), cigar_mismatch (one base),
 *   cigar_ins (bases only in the reference: SAM 'D') and cigar_del (bases
 *   only in the read: SAM 'I', or 'S' with cigar_attr_trim).
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     matches_cigar_vector_append_{insertion,deletion,match,mismatch,cigar_element}(vector_t*,uint64_t*,...)
 *     matches_cigar_buffer_add_cigar_element[__attr](cigar_element_t**,cigar_t,uint64_t[,cigar_attr_t])
 *     matches_cigar_buffer_add_mismatch(cigar_element_t**,uint8_t)
 *     matches_cigar_vector_insert_cigar_element(vector_t*,int64_t,cigar_element_t*)
 *     matches_cigar_compute_{event_distance,edit_distance,edit_distance__excluding_clipping,
 *         matching_bases,effective_length}(vector_t*,uint64_t offset,uint64_t length)
 *     matches_cigar_element_effective_length(const cigar_element_t*)
 *     matches_cigar_compute_error_quality(vector_t*,uint64_t,uint64_t,uint8_t* quality_mask,uint64_t)
 *     match_cigar_print(FILE*,vector_t*,uint64_t,uint64_t)
 *     matches_counters_{new,clear,delete,add,sub,get_count,get_total_count,get_num_counters,
 *         get_counts,compact,compute_matches_to_report,count_first_subdominant}
 *     matches_metrics_{init,update,set_proper_length,set_read_length,set_swg_match_score,
 *         set_region_profile_metrics,add_accepted_candidates}
 *     matches_{new,clear,delete,is_mapped,recompute_metrics,update_mcs,update_limited_exact_matches,
 *         index_match,lookup_match,get_num_match_traces,get_num_match_traces_extended,
 *         get_primary_match,get_match_traces,traces_sort_by_genomic_position,
 *         add_match_trace_insert_sorted,add_match_trace_preserve_sorted,
 *         add_match_trace_local_pending,local_pending_add_to_extended_matches}
 *     match_trace_copy(match_trace_t* dst,match_trace_t* src) (no prototype in a header)
 *     match_trace_{get_cigar_buffer,get_cigar_length,cmp_swg_score,cmp_genomic_position}
 *     matches_classify, matches_classify_{pattern_viable,max_matches_searched,min_depth_searched,
 *         neighbourhood_fallback,local_alignment_fallback,adjust_max_error_by_strata_after_best,
 *         compute_max_search_error}
 *     matches_predictors_{compute_unique,compute_mmaps,compute_se}
 *     select_parameters_{init,instantiate_values}, archive_select_configure_{se,pe}(archive_search_t*)
 *     archive_select_se_matches(select_parameters_t*,matches_t*), archive_select_se_matches_discard
 *     archive_score_probability_to_mapq(double), archive_score_logit, archive_score_logit_{unique,mmap}
 *     archive_score_matches_se_get_model, archive_score_matches_se_logit_{unique,mmap}
 *     archive_score_matches_se[_default,_classify](archive_search_t*,matches_t*)
 *   and builds or reads these structs field by field: cigar_element_t;
 *   match_trace_t (type, scores and distances, error_quality, sequence_name,
 *   text_position, strand, bs_strand, mapq_score, match_alignment); matches_t
 *   (counters, metrics, classification, max_complete_stratum,
 *   limited_exact_matches, cigar_vector, match_traces_local,
 *   match_traces_extended); matches_metrics_t; matches_predictors_t;
 *   select_parameters_t; search_parameters_t (mapping_mode,
 *   select_parameters, match_alignment_model, mapq_model_se);
 *   archive_search_t.search_parameters; archive_score_logit_model_t.
 */

#include "gem3_test.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"
#include "matches/matches_counters.h"
#include "matches/match_trace.h"
#include "matches/classify/matches_classify.h"
#include "matches/classify/matches_metrics.h"
#include "matches/classify/matches_predictors.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_select.h"
#include "archive/search/archive_select_parameters.h"
#include "archive/score/archive_score_se.h"
#include "archive/score/archive_score_logit.h"
#include "archive/score/archive_score_logit_models.h"
#include <math.h>

/* Defined in src/ without a prototype in a header */
void match_trace_copy(match_trace_t* const match_trace_dst,match_trace_t* const match_trace_src);
void matches_add_match_trace_insert_sorted(matches_t* const matches,match_trace_t* const match_trace);
void matches_add_match_trace_preserve_sorted(matches_t* const matches);
void matches_index_match(matches_t* const matches,match_trace_t* const match_trace);
match_trace_t* matches_lookup_match(
    matches_t* const matches,const uint64_t begin_position,const uint64_t effective_length);
void matches_predictors_compute_unique(
    matches_predictors_t* const predictors,matches_metrics_t* const metrics,
    const uint64_t primary_edit_distance,const uint64_t primary_event_distance,
    const int32_t primary_swg_score,const double primary_error_quality,
    const double primary_template_size_sigma);
void matches_predictors_compute_mmaps(
    matches_predictors_t* const predictors,matches_metrics_t* const metrics,
    const uint64_t primary_edit_distance,const uint64_t primary_event_distance,
    const int32_t primary_swg_score,const double primary_error_quality,
    const double primary_template_size_sigma);
void archive_select_se_matches_discard(matches_t* const matches,const uint64_t strata_to_report);
const archive_score_logit_model_t* archive_score_matches_se_get_model(search_parameters_t* const search_parameters);

#define MAX_ELEMENTS 512

/*
 * CIGAR reference model: a list of (type,length,attributes); a mismatch has
 * length 1 and keeps its base in 'mismatch'
 */
typedef struct {
  cigar_t type;
  int64_t length;
  cigar_attr_t attributes;
  uint8_t mismatch;
} ref_elm_t;

static void random_elements(gt_rng_t* const rng,ref_elm_t* const elms,const uint64_t n) {
  uint64_t i;
  for (i=0;i<n;++i) {
    const uint64_t r = gt_rng_below(rng,10);
    elms[i].type = (r<4) ? cigar_match : (r<7) ? cigar_mismatch : (r<8) ? cigar_ins : cigar_del;
    elms[i].length = (elms[i].type==cigar_mismatch) ? 1 : (int64_t)gt_rng_range(rng,1,30);
    elms[i].attributes = (cigar_attr_t)gt_rng_below(rng,3);
    elms[i].mismatch = (uint8_t)gt_rng_below(rng,4);
  }
}
// Run-length merge as the appenders describe it: consecutive match, ins or
// del elements merge (the merged element loses its attributes); mismatches
// never merge
static uint64_t ref_merge(const ref_elm_t* const in,const uint64_t n,ref_elm_t* const out) {
  uint64_t i, m = 0;
  for (i=0;i<n;++i) {
    if (m>0 && in[i].type!=cigar_mismatch && out[m-1].type==in[i].type) {
      out[m-1].length += in[i].length;
      out[m-1].attributes = cigar_attr_none;
    } else {
      out[m++] = in[i];
    }
  }
  return m;
}
static void to_cigar_element(const ref_elm_t* const e,cigar_element_t* const c) {
  memset(c,0,sizeof(*c));
  c->type = e->type;
  c->attributes = e->attributes;
  if (e->type==cigar_mismatch) c->mismatch = e->mismatch;
  else c->length = (int32_t)e->length;
}
static bool same_element(const cigar_element_t* const c,const ref_elm_t* const e) {
  if (c->type!=e->type || c->attributes!=e->attributes) return false;
  if (e->type==cigar_mismatch) return c->mismatch==e->mismatch;
  return c->length==e->length;
}
// Load reference elements into a vector (no merging) at its end
static void load_vector(vector_t* const v,const ref_elm_t* const elms,const uint64_t n) {
  uint64_t i;
  vector_reserve_additional(v,n);
  for (i=0;i<n;++i) to_cigar_element(elms+i,vector_get_mem(v,cigar_element_t)+vector_get_used(v)+i);
  vector_add_used(v,n);
}

/*
 * Appending merges runs of the same operation
 */
GT_TEST(cigar_append_merges_runs) {
  gt_rng_t rng = gt_rng_new(101);
  static ref_elm_t in[MAX_ELEMENTS], expected[MAX_ELEMENTS];
  uint64_t trial;
  for (trial=0;trial<300;++trial) {
    const uint64_t n = gt_rng_range(&rng,1,60);
    random_elements(&rng,in,n);
    const uint64_t m = ref_merge(in,n,expected);
    vector_t* const v = vector_new(4,cigar_element_t);
    uint64_t cigar_length = 0, i;
    for (i=0;i<n;++i) {
      cigar_element_t c;
      to_cigar_element(in+i,&c);
      switch (trial%2) {
        case 0: matches_cigar_vector_append_cigar_element(v,&cigar_length,&c); break;
        default: // The typed appenders
          switch (c.type) {
            case cigar_match: matches_cigar_vector_append_match(v,&cigar_length,c.length,c.attributes); break;
            case cigar_mismatch: matches_cigar_vector_append_mismatch(v,&cigar_length,c.mismatch,c.attributes); break;
            case cigar_ins: matches_cigar_vector_append_insertion(v,&cigar_length,c.length,c.attributes); break;
            default: matches_cigar_vector_append_deletion(v,&cigar_length,c.length,c.attributes); break;
          }
          break;
      }
    }
    GT_CHECK_EQ_UINT(cigar_length,m);
    GT_CHECK_EQ_UINT(vector_get_used(v),m);
    for (i=0;i<m && i<vector_get_used(v);++i) {
      GT_CHECK_MSG(same_element(vector_get_mem(v,cigar_element_t)+i,expected+i),
          "trial %" PRIu64 " element %" PRIu64 " differs",trial,i);
    }
    vector_delete(v);
  }
}

/*
 * A new CIGAR (length 0) never merges into the previous CIGAR stored in the
 * same vector; a CIGAR that has begun does
 */
GT_TEST(cigar_append_starts_new_cigar) {
  vector_t* const v = vector_new(4,cigar_element_t);
  uint64_t first = 0, second = 0;
  matches_cigar_vector_append_match(v,&first,10,cigar_attr_none);
  matches_cigar_vector_append_insertion(v,&first,2,cigar_attr_none);
  matches_cigar_vector_append_match(v,&first,5,cigar_attr_none);
  GT_CHECK_EQ_UINT(first,3);
  // Second CIGAR starts with a match: a new element, not 5+7
  matches_cigar_vector_append_match(v,&second,7,cigar_attr_homopolymer);
  GT_CHECK_EQ_UINT(second,1);
  GT_REQUIRE(vector_get_used(v)==4);
  cigar_element_t* const e = vector_get_mem(v,cigar_element_t);
  GT_CHECK_EQ_INT(e[2].length,5);
  GT_CHECK_EQ_INT(e[3].type,cigar_match);
  GT_CHECK_EQ_INT(e[3].length,7);
  GT_CHECK_EQ_INT(e[3].attributes,cigar_attr_homopolymer);
  // ... and then merges within itself, dropping the attribute
  matches_cigar_vector_append_match(v,&second,3,cigar_attr_trim);
  GT_CHECK_EQ_UINT(second,1);
  GT_CHECK_EQ_INT(vector_get_mem(v,cigar_element_t)[3].length,10);
  GT_CHECK_EQ_INT(vector_get_mem(v,cigar_element_t)[3].attributes,cigar_attr_none);
  // Same for insertions and deletions
  uint64_t third = 0;
  matches_cigar_vector_append_deletion(v,&third,4,cigar_attr_trim);
  matches_cigar_vector_append_deletion(v,&third,1,cigar_attr_trim);
  GT_CHECK_EQ_UINT(third,1);
  GT_CHECK_EQ_INT(vector_get_last_elm(v,cigar_element_t)->length,5);
  GT_CHECK_EQ_INT(vector_get_last_elm(v,cigar_element_t)->attributes,cigar_attr_none);
  // Mismatches never merge
  matches_cigar_vector_append_mismatch(v,&third,1,cigar_attr_none);
  matches_cigar_vector_append_mismatch(v,&third,1,cigar_attr_none);
  GT_CHECK_EQ_UINT(third,3);
  vector_delete(v);
}

/*
 * The buffer appender works through a sentinel pointing at the last element
 * written (type cigar_null before the first one)
 */
GT_TEST(cigar_buffer_sentinel) {
  gt_rng_t rng = gt_rng_new(102);
  static ref_elm_t in[MAX_ELEMENTS], expected[MAX_ELEMENTS];
  static cigar_element_t buffer[MAX_ELEMENTS+1];
  uint64_t trial;
  for (trial=0;trial<300;++trial) {
    const uint64_t n = gt_rng_range(&rng,1,60);
    random_elements(&rng,in,n);
    const uint64_t m = ref_merge(in,n,expected);
    memset(buffer,0,sizeof(buffer));
    buffer[0].type = cigar_null;
    cigar_element_t* sentinel = buffer;
    uint64_t i;
    for (i=0;i<n;++i) {
      if (in[i].type==cigar_mismatch) {
        matches_cigar_buffer_add_mismatch(&sentinel,in[i].mismatch);
      } else if (in[i].attributes==cigar_attr_none && (trial%2)) {
        matches_cigar_buffer_add_cigar_element(&sentinel,in[i].type,in[i].length);
      } else {
        matches_cigar_buffer_add_cigar_element__attr(&sentinel,in[i].type,in[i].length,in[i].attributes);
      }
    }
    GT_CHECK_EQ_UINT((uint64_t)(sentinel-buffer)+1,m);
    for (i=0;i<m;++i) {
      ref_elm_t e = expected[i];
      if (e.type==cigar_mismatch) e.attributes = cigar_attr_none; // add_mismatch stores no attributes
      GT_CHECK_MSG(same_element(buffer+i,&e),"trial %" PRIu64 " element %" PRIu64 " differs",trial,i);
    }
  }
}

/*
 * Inserting an element at a position shifts the rest right
 */
GT_TEST(cigar_insert_element) {
  gt_rng_t rng = gt_rng_new(103);
  static ref_elm_t in[MAX_ELEMENTS];
  uint64_t trial;
  for (trial=0;trial<200;++trial) {
    const uint64_t n = gt_rng_range(&rng,0,40);
    random_elements(&rng,in,n+1);
    vector_t* const v = vector_new(2,cigar_element_t);
    load_vector(v,in,n);
    const uint64_t pos = gt_rng_below(&rng,n+1);
    cigar_element_t c;
    to_cigar_element(in+n,&c);
    matches_cigar_vector_insert_cigar_element(v,pos,&c);
    GT_REQUIRE(vector_get_used(v)==n+1);
    uint64_t i;
    for (i=0;i<=n;++i) {
      const ref_elm_t* const e = (i<pos) ? in+i : (i==pos) ? in+n : in+i-1;
      GT_CHECK_MSG(same_element(vector_get_mem(v,cigar_element_t)+i,e),
          "trial %" PRIu64 ": element %" PRIu64 " after inserting at %" PRIu64,trial,i,pos);
    }
    vector_delete(v);
  }
}

/*
 * Distances and lengths of random CIGARs (stored after an unrelated prefix,
 * so offsets matter), recomputed by expanding each CIGAR into columns
 */
GT_TEST(cigar_distances) {
  gt_rng_t rng = gt_rng_new(104);
  static ref_elm_t prefix[MAX_ELEMENTS], in[MAX_ELEMENTS];
  uint64_t trial;
  for (trial=0;trial<500;++trial) {
    const uint64_t np = gt_rng_below(&rng,8), n = gt_rng_range(&rng,1,40);
    random_elements(&rng,prefix,np);
    random_elements(&rng,in,n);
    vector_t* const v = vector_new(2,cigar_element_t);
    load_vector(v,prefix,np);
    load_vector(v,in,n);
    uint64_t events = 0, edits = 0, edits_noclip = 0, matching = 0, ref_span = 0, i, k;
    for (i=0;i<n;++i) {
      // Column by column: a match column matches, a mismatch column is one
      // edit, each indel column is one edit, and an indel run is one event
      for (k=0;k<(uint64_t)in[i].length;++k) {
        switch (in[i].type) {
          case cigar_match: ++matching; ++ref_span; break;
          case cigar_mismatch: ++edits; ++edits_noclip; ++ref_span; break;
          case cigar_ins: ++edits; ++edits_noclip; ++ref_span; break; // reference-only base
          default: ++edits; if (i>0 && i<n-1) ++edits_noclip; break; // read-only base; at an end it is a trim
        }
      }
      if (in[i].type!=cigar_match) ++events;
      GT_CHECK_EQ_INT(matches_cigar_element_effective_length(vector_get_mem(v,cigar_element_t)+np+i),
          (in[i].type==cigar_del) ? 0 : in[i].length);
    }
    GT_CHECK_EQ_UINT(matches_cigar_compute_event_distance(v,np,n),events);
    GT_CHECK_EQ_UINT(matches_cigar_compute_edit_distance(v,np,n),edits);
    GT_CHECK_EQ_UINT(matches_cigar_compute_edit_distance__excluding_clipping(v,np,n),edits_noclip);
    GT_CHECK_EQ_UINT(matches_cigar_compute_matching_bases(v,np,n),matching);
    GT_CHECK_EQ_INT(matches_cigar_compute_effective_length(v,np,n),ref_span);
    vector_delete(v);
  }
  // Hand-worked: 10M 1X 3(ref-only) 5M 2(read-only) = 1+3+2 edits, 3 events, span 19
  ref_elm_t hand[] = {
    {cigar_match,10,cigar_attr_none,0}, {cigar_mismatch,1,cigar_attr_none,2},
    {cigar_ins,3,cigar_attr_none,0}, {cigar_match,5,cigar_attr_none,0},
    {cigar_del,2,cigar_attr_none,0},
  };
  vector_t* const v = vector_new(2,cigar_element_t);
  load_vector(v,hand,5);
  GT_CHECK_EQ_UINT(matches_cigar_compute_edit_distance(v,0,5),6);
  GT_CHECK_EQ_UINT(matches_cigar_compute_edit_distance__excluding_clipping(v,0,5),4);
  GT_CHECK_EQ_UINT(matches_cigar_compute_event_distance(v,0,5),3);
  GT_CHECK_EQ_UINT(matches_cigar_compute_matching_bases(v,0,5),15);
  GT_CHECK_EQ_INT(matches_cigar_compute_effective_length(v,0,5),19);
  vector_delete(v);
}

/*
 * Average quality of the error columns: a mismatch takes the quality of its
 * read base, a reference-only run (inside the alignment) the quality of the
 * read base that follows it once, a read-only run (inside) each of its
 * bases; runs at either end are not errors. No mask or no CIGAR gives
 * SEQUENCE_QUALITIES_MAX.
 */
static float ref_error_quality(const ref_elm_t* const in,const uint64_t n,const uint8_t* const q) {
  uint64_t i, k, read_pos = 0, count = 0, sum = 0;
  for (i=0;i<n;++i) {
    const bool inside = (i>0 && i<n-1);
    switch (in[i].type) {
      case cigar_match: read_pos += in[i].length; break;
      case cigar_mismatch: sum += q[read_pos++]; ++count; break;
      case cigar_ins: if (inside) { sum += q[read_pos]; ++count; } break;
      default:
        for (k=0;k<(uint64_t)in[i].length;++k,++read_pos) {
          if (inside) { sum += q[read_pos]; ++count; }
        }
        break;
    }
  }
  return count ? (float)sum/(float)count : (float)SEQUENCE_QUALITIES_MAX;
}
GT_TEST(cigar_error_quality) {
  gt_rng_t rng = gt_rng_new(105);
  static ref_elm_t in[MAX_ELEMENTS];
  static uint8_t q[4096];
  uint64_t trial, i;
  for (trial=0;trial<500;++trial) {
    const uint64_t n = gt_rng_range(&rng,1,30);
    random_elements(&rng,in,n);
    uint64_t read_length = 0;
    for (i=0;i<n;++i) if (in[i].type!=cigar_ins) read_length += in[i].length;
    for (i=0;i<read_length+1;++i) q[i] = (uint8_t)gt_rng_range(&rng,2,41);
    vector_t* const v = vector_new(2,cigar_element_t);
    load_vector(v,in,n);
    const float got = matches_cigar_compute_error_quality(v,0,n,q,read_length);
    const float want = ref_error_quality(in,n,q);
    GT_CHECK_MSG(fabsf(got-want) < 1e-4f,"trial %" PRIu64 ": %f != %f",trial,got,want);
    GT_CHECK_EQ_INT((int)matches_cigar_compute_error_quality(v,0,n,NULL,0),SEQUENCE_QUALITIES_MAX);
    GT_CHECK_EQ_INT((int)matches_cigar_compute_error_quality(v,0,0,q,read_length),SEQUENCE_QUALITIES_MAX);
    vector_delete(v);
  }
  // Hand-worked: 2M 1X(q=30) 3M 1X(q=10) -> (30+10)/2
  ref_elm_t hand[] = {
    {cigar_match,2,cigar_attr_none,0}, {cigar_mismatch,1,cigar_attr_none,0},
    {cigar_match,3,cigar_attr_none,0}, {cigar_mismatch,1,cigar_attr_none,0},
  };
  uint8_t hq[7] = {40,40,30,40,40,40,10};
  vector_t* const v = vector_new(2,cigar_element_t);
  load_vector(v,hand,4);
  GT_CHECK(fabsf(matches_cigar_compute_error_quality(v,0,4,hq,7)-20.0f) < 1e-6f);
  vector_delete(v);
}

/*
 * The debug printer: M, 1X, I/i, D/d/S
 */
GT_TEST(cigar_print) {
  ref_elm_t hand[] = {
    {cigar_match,12,cigar_attr_none,0}, {cigar_mismatch,1,cigar_attr_none,1},
    {cigar_ins,2,cigar_attr_none,0}, {cigar_ins,3,cigar_attr_homopolymer,0},
    {cigar_del,4,cigar_attr_none,0}, {cigar_del,5,cigar_attr_homopolymer,0},
    {cigar_del,6,cigar_attr_trim,0}, {cigar_match,7,cigar_attr_trim,0},
  };
  vector_t* const v = vector_new(2,cigar_element_t);
  load_vector(v,hand,1); // offset 1 below
  load_vector(v,hand,8);
  char buffer[256];
  FILE* const stream = fmemopen(buffer,sizeof(buffer),"w");
  GT_REQUIRE(stream!=NULL);
  match_cigar_print(stream,v,1,8);
  fclose(stream);
  GT_CHECK_EQ_STR(buffer,"12M1X2I3i4D5d6S7M");
  vector_delete(v);
}

/*
 * Counters: one per edit distance, with a running total
 */
GT_TEST(counters_random_ops) {
  gt_rng_t rng = gt_rng_new(106);
  uint64_t ref[64];
  uint64_t trial;
  for (trial=0;trial<200;++trial) {
    matches_counters_t* const c = matches_counters_new();
    memset(ref,0,sizeof(ref));
    uint64_t used = 0, total = 0, op, d;
    for (op=0;op<50;++op) {
      d = gt_rng_below(&rng,40);
      if (gt_rng_chance(&rng,1,4) && d<used && ref[d]>0) {
        matches_counters_sub(c,d,1);
        --ref[d]; --total;
      } else {
        const uint64_t k = gt_rng_range(&rng,1,3);
        matches_counters_add(c,d,k);
        ref[d] += k; total += k;
        if (d+1>used) used = d+1;
      }
    }
    GT_CHECK_EQ_UINT(matches_counters_get_num_counters(c),used);
    GT_CHECK_EQ_UINT(matches_counters_get_total_count(c),total);
    for (d=0;d<used;++d) {
      GT_CHECK_EQ_UINT(matches_counters_get_count(c,d),ref[d]);
      GT_CHECK_EQ_UINT(matches_counters_get_counts(c)[d],ref[d]);
    }
    // First subdominant: the count of the first non-empty stratum after the best one
    uint64_t best = 0, sub = 0;
    while (best<used && ref[best]==0) ++best;
    for (d=best+1;d<used;++d) if (ref[d]) { sub = ref[d]; break; }
    GT_CHECK_EQ_UINT(matches_counters_count_first_subdominant(c),sub);
    // Compact drops trailing empty strata
    uint64_t last = used;
    while (last>0 && ref[last-1]==0) --last;
    GT_CHECK_EQ_UINT(matches_counters_compact(c),last);
    GT_CHECK_EQ_UINT(matches_counters_get_num_counters(c),last);
    matches_counters_clear(c);
    GT_CHECK_EQ_UINT(matches_counters_get_num_counters(c),0);
    GT_CHECK_EQ_UINT(matches_counters_get_total_count(c),0);
    GT_CHECK_EQ_UINT(matches_counters_count_first_subdominant(c),0);
    matches_counters_delete(c);
  }
}

/*
 * Matches to report for -m (strata) and -M (matches), hand-worked:
 * at least min_reported_strata strata are taken, counting from the first
 * non-empty one; further strata are added while the total stays within
 * max_reported_matches. strata_to_report counts the strata taken from 0.
 */
GT_TEST(counters_matches_to_report) {
  static const struct {
    uint64_t counts[8];
    uint64_t num_counts;
    uint64_t min_strata, max_matches;
    uint64_t matches, strata;
  } cases[] = {
    { {0},            0, 1, 5, 0, 0 },  // nothing found
    { {0,0,0},        3, 1, 5, 0, 0 },  // only empty strata
    { {1,2,3,4},      4, 1, 5, 3, 2 },  // 0:1, then 1:2 (3 <= 5); 2:3 would make 6
    { {1,2,3,4},      4, 2, 5, 3, 2 },  // two strata are 0 and 1
    { {1,2,3,4},      4, 3, 5, 6, 3 },  // three strata even beyond -M
    { {1,2,3,4},      4, 1, 100, 10, 4 }, // everything fits
    { {0,0,6,1},      4, 1, 5, 6, 3 },  // first non-empty stratum (2) is complete even past -M
    { {0,2,0,0,3},    5, 1, 5, 5, 5 },  // empty strata between are taken as they are free
    { {0,2,0,0,3},    5, 1, 4, 2, 4 },  // 2 + 3 > 4: stops before stratum 4
    { {7},            1, 1, 5, 7, 1 },
  };
  uint64_t i, j;
  for (i=0;i<sizeof(cases)/sizeof(cases[0]);++i) {
    matches_counters_t* const c = matches_counters_new();
    for (j=0;j<cases[i].num_counts;++j) {
      if (cases[i].counts[j]) matches_counters_add(c,j,cases[i].counts[j]);
    }
    uint64_t matches = 99, strata = 99;
    matches_counters_compute_matches_to_report(c,cases[i].min_strata,cases[i].max_matches,&matches,&strata);
    GT_CHECK_MSG(matches==cases[i].matches && strata==cases[i].strata,
        "case %" PRIu64 ": (%" PRIu64 ",%" PRIu64 ") expected (%" PRIu64 ",%" PRIu64 ")",
        i,matches,strata,cases[i].matches,cases[i].strata);
    matches_counters_delete(c);
  }
}

/*
 * Metrics keep the two smallest event and edit distances and the two
 * largest SWG scores of the multiset of updates
 */
static int cmp_u64(const void* a,const void* b) {
  const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
  return (x<y) ? -1 : (x>y);
}
static int cmp_i32_desc(const void* a,const void* b) {
  const int32_t x = *(const int32_t*)a, y = *(const int32_t*)b;
  return (x>y) ? -1 : (x<y);
}
GT_TEST(metrics_update) {
  gt_rng_t rng = gt_rng_new(107);
  uint64_t trial;
  for (trial=0;trial<500;++trial) {
    matches_metrics_t m;
    matches_metrics_init(&m);
    GT_CHECK_EQ_UINT(m.min1_edit_distance,UINT32_MAX);
    GT_CHECK_EQ_UINT(m.min2_event_distance,UINT32_MAX);
    GT_CHECK_EQ_INT(m.max1_swg_score,INT32_MIN);
    GT_CHECK_EQ_UINT(m.candidates_accepted,0);
    GT_CHECK_EQ_UINT(m.mapq,0);
    const uint64_t n = gt_rng_range(&rng,1,12);
    uint64_t ev[16], ed[16], i;
    int32_t sw[16];
    for (i=0;i<n;++i) {
      ev[i] = gt_rng_below(&rng,6);
      ed[i] = gt_rng_below(&rng,6);
      sw[i] = (int32_t)gt_rng_below(&rng,8) - 2;
      matches_metrics_update(&m,ev[i],ed[i],sw[i]);
    }
    qsort(ev,n,sizeof(uint64_t),cmp_u64);
    qsort(ed,n,sizeof(uint64_t),cmp_u64);
    qsort(sw,n,sizeof(int32_t),cmp_i32_desc);
    GT_CHECK_EQ_UINT(m.min1_event_distance,ev[0]);
    GT_CHECK_EQ_UINT(m.min1_edit_distance,ed[0]);
    GT_CHECK_EQ_INT(m.max1_swg_score,sw[0]);
    GT_CHECK_EQ_UINT(m.min2_event_distance,n>1 ? ev[1] : UINT32_MAX);
    GT_CHECK_EQ_UINT(m.min2_edit_distance,n>1 ? ed[1] : UINT32_MAX);
    GT_CHECK_EQ_INT(m.max2_swg_score,n>1 ? sw[1] : INT32_MIN);
  }
  matches_metrics_t m;
  matches_metrics_init(&m);
  matches_metrics_set_proper_length(&m,9.5);
  matches_metrics_set_read_length(&m,150);
  matches_metrics_set_swg_match_score(&m,2);
  matches_metrics_set_region_profile_metrics(&m,20,31,0.25);
  matches_metrics_add_accepted_candidates(&m,3);
  matches_metrics_add_accepted_candidates(&m,4);
  GT_CHECK(m.proper_length==9.5);
  GT_CHECK_EQ_UINT(m.read_length,150);
  GT_CHECK_EQ_INT(m.swg_match_score,2);
  GT_CHECK_EQ_UINT(m.avg_region_length,20);
  GT_CHECK_EQ_UINT(m.max_region_length,31);
  GT_CHECK(m.kmer_frequency==0.25);
  GT_CHECK_EQ_UINT(m.candidates_accepted,7);
}

/*
 * Match traces built by hand (no locator needed)
 */
static char seq_names[3][8] = {"chrA","chrB","chrC"};
static void make_trace(
    match_trace_t* const t,const int32_t swg,const uint64_t event,const uint64_t edit,
    const float quality,const int name,const uint64_t position) {
  memset(t,0,sizeof(*t));
  t->type = match_type_regular;
  t->swg_score = swg;
  t->event_distance = event;
  t->edit_distance = edit;
  t->error_quality = quality;
  t->sequence_name = seq_names[name];
  t->text_position = position;
  t->strand = Forward;
  t->match_alignment.match_position = position*2+(uint64_t)name*1000000;
  t->match_alignment.effective_length = 100;
}
// The documented order: SWG score high to low, then event distance, edit
// distance low to high, error quality high to low, then position
static int ref_cmp_swg(const match_trace_t* const a,const match_trace_t* const b) {
  if (a->swg_score != b->swg_score) return (a->swg_score > b->swg_score) ? -1 : 1;
  if (a->event_distance != b->event_distance) return (a->event_distance < b->event_distance) ? -1 : 1;
  if (a->edit_distance != b->edit_distance) return (a->edit_distance < b->edit_distance) ? -1 : 1;
  if (a->error_quality != b->error_quality) return (a->error_quality > b->error_quality) ? -1 : 1;
  if (a->sequence_name != b->sequence_name) return (a->sequence_name < b->sequence_name) ? -1 : 1;
  if (a->text_position != b->text_position) return (a->text_position < b->text_position) ? -1 : 1;
  return 0;
}
static int sign(const int x) { return (x>0)-(x<0); }
static void random_trace(gt_rng_t* const rng,match_trace_t* const t) {
  make_trace(t,(int32_t)gt_rng_range(rng,60,100),gt_rng_below(rng,5),gt_rng_below(rng,5),
      (float)gt_rng_range(rng,20,40),(int)gt_rng_below(rng,3),gt_rng_below(rng,5000));
}

GT_TEST(matches_insert_sorted) {
  gt_rng_t rng = gt_rng_new(108);
  static match_trace_t traces[64];
  uint64_t trial, i;
  for (trial=0;trial<200;++trial) {
    matches_t* const matches = matches_new();
    GT_CHECK(!matches_is_mapped(matches));
    const uint64_t n = gt_rng_range(&rng,1,40);
    for (i=0;i<n;++i) {
      random_trace(&rng,traces+i);
      matches_add_match_trace_insert_sorted(matches,traces+i);
    }
    GT_REQUIRE(matches_get_num_match_traces(matches)==n);
    GT_CHECK(matches_is_mapped(matches));
    match_trace_t** const mt = matches_get_match_traces(matches);
    for (i=0;i+1<n;++i) {
      GT_CHECK_MSG(ref_cmp_swg(mt[i],mt[i+1]) <= 0,"trial %" PRIu64 ": not sorted at %" PRIu64,trial,i);
      // Equal keys keep insertion order
      if (ref_cmp_swg(mt[i],mt[i+1])==0) GT_CHECK(mt[i] < mt[i+1]);
    }
    for (i=0;i<n;++i) {
      const uint64_t j = gt_rng_below(&rng,n);
      GT_CHECK_EQ_INT(sign(match_trace_cmp_swg_score((const match_trace_t**)&mt[i],(const match_trace_t**)&mt[j])),
          sign(ref_cmp_swg(mt[i],mt[j])));
    }
    GT_CHECK(matches_get_primary_match(matches)==mt[0]);
    matches_delete(matches);
  }
}

GT_TEST(matches_preserve_sorted) {
  gt_rng_t rng = gt_rng_new(109);
  static match_trace_t traces[64];
  uint64_t trial, i;
  for (trial=0;trial<200;++trial) {
    matches_t* const matches = matches_new();
    const uint64_t n = gt_rng_range(&rng,1,30);
    for (i=0;i<n;++i) {
      random_trace(&rng,traces+i);
      matches_add_match_trace_insert_sorted(matches,traces+i);
    }
    // Improve or worsen one match in place, then restore the order
    match_trace_t** const mt = matches_get_match_traces(matches);
    match_trace_t* const changed = mt[gt_rng_below(&rng,n)];
    changed->swg_score = (int32_t)gt_rng_range(&rng,50,110);
    matches_add_match_trace_preserve_sorted(matches);
    for (i=0;i+1<n;++i) {
      GT_CHECK_MSG(ref_cmp_swg(mt[i],mt[i+1]) <= 0,"trial %" PRIu64 ": not sorted at %" PRIu64,trial,i);
    }
    uint64_t found = 0;
    for (i=0;i<n;++i) found += (mt[i]==changed);
    GT_CHECK_EQ_UINT(found,1);
    matches_delete(matches);
  }
}

GT_TEST(matches_sort_by_genomic_position) {
  gt_rng_t rng = gt_rng_new(110);
  static match_trace_t traces[64];
  match_trace_t* ptrs[64];
  uint64_t trial, i;
  for (trial=0;trial<200;++trial) {
    const uint64_t n = gt_rng_range(&rng,1,50);
    for (i=0;i<n;++i) {
      // Distinct (name,position) keys; names compared as strings, not pointers
      make_trace(traces+i,80,0,0,30.0f,(int)gt_rng_below(&rng,3),gt_rng_below(&rng,100)*64+i);
      ptrs[i] = traces+i;
    }
    matches_traces_sort_by_genomic_position(ptrs,n);
    for (i=0;i+1<n;++i) {
      const int c = strcmp(ptrs[i]->sequence_name,ptrs[i+1]->sequence_name);
      GT_CHECK(c<0 || (c==0 && ptrs[i]->text_position < ptrs[i+1]->text_position));
      GT_CHECK(match_trace_cmp_genomic_position((const match_trace_t**)&ptrs[i],(const match_trace_t**)&ptrs[i+1]) < 0);
      GT_CHECK(match_trace_cmp_genomic_position((const match_trace_t**)&ptrs[i+1],(const match_trace_t**)&ptrs[i]) > 0);
    }
    GT_CHECK(match_trace_cmp_genomic_position((const match_trace_t**)&ptrs[0],(const match_trace_t**)&ptrs[0])==0);
  }
}

/*
 * Index by begin and by end (begin+effective length) positions
 */
GT_TEST(matches_index_lookup) {
  matches_t* const matches = matches_new();
  static match_trace_t t[3];
  make_trace(t+0,90,1,1,30.0f,0,100);
  make_trace(t+1,80,2,2,30.0f,0,500);
  make_trace(t+2,70,3,3,30.0f,1,900);
  t[0].match_alignment.match_position = 1000; t[0].match_alignment.effective_length = 100;
  t[1].match_alignment.match_position = 5000; t[1].match_alignment.effective_length = 98;
  t[2].match_alignment.match_position = 9000; t[2].match_alignment.effective_length = 103;
  uint64_t i;
  for (i=0;i<3;++i) matches_index_match(matches,t+i);
  GT_CHECK(matches_lookup_match(matches,1000,100)==t+0);
  GT_CHECK(matches_lookup_match(matches,1000,7)==t+0);     // same begin
  GT_CHECK(matches_lookup_match(matches,4998,100)==t+1);   // same end (5098)
  GT_CHECK(matches_lookup_match(matches,9003,100)==t+2);   // same end (9103)
  GT_CHECK(matches_lookup_match(matches,9001,100)==NULL);
  GT_CHECK(matches_lookup_match(matches,0,0)==NULL);
  matches_clear(matches);
  GT_CHECK(matches_lookup_match(matches,1000,100)==NULL);
  matches_delete(matches);
}

GT_TEST(matches_strata_bookkeeping) {
  matches_t* const matches = matches_new();
  GT_CHECK_EQ_UINT(matches->max_complete_stratum,0);
  matches_update_mcs(matches,3);
  matches_update_mcs(matches,2);
  GT_CHECK_EQ_UINT(matches->max_complete_stratum,3);
  matches_update_limited_exact_matches(matches,0);
  GT_CHECK(!matches->limited_exact_matches);
  GT_CHECK_EQ_UINT(matches->max_complete_stratum,3);
  matches_update_limited_exact_matches(matches,4);
  GT_CHECK(matches->limited_exact_matches);
  GT_CHECK_EQ_UINT(matches->max_complete_stratum,0);
  GT_CHECK_EQ_UINT(matches_counters_get_count(matches->counters,0),4);
  matches_clear(matches);
  GT_CHECK(!matches->limited_exact_matches);
  GT_CHECK_EQ_UINT(matches_counters_get_total_count(matches->counters),0);
  matches_delete(matches);
}

GT_TEST(matches_recompute_and_accessors) {
  gt_rng_t rng = gt_rng_new(111);
  static match_trace_t traces[32];
  matches_t* const matches = matches_new();
  matches_metrics_t expected;
  matches_metrics_init(&expected);
  uint64_t i;
  for (i=0;i<20;++i) {
    random_trace(&rng,traces+i);
    matches_add_match_trace_insert_sorted(matches,traces+i);
    matches_metrics_update(&expected,traces[i].event_distance,traces[i].edit_distance,traces[i].swg_score);
  }
  matches_recompute_metrics(matches);
  GT_CHECK_EQ_UINT(matches->metrics.min1_edit_distance,expected.min1_edit_distance);
  GT_CHECK_EQ_UINT(matches->metrics.min2_edit_distance,expected.min2_edit_distance);
  GT_CHECK_EQ_UINT(matches->metrics.min1_event_distance,expected.min1_event_distance);
  GT_CHECK_EQ_UINT(matches->metrics.min2_event_distance,expected.min2_event_distance);
  GT_CHECK_EQ_INT(matches->metrics.max1_swg_score,expected.max1_swg_score);
  GT_CHECK_EQ_INT(matches->metrics.max2_swg_score,expected.max2_swg_score);
  // CIGAR accessors index the matches' CIGAR vector
  ref_elm_t hand[] = { {cigar_match,40,cigar_attr_none,0}, {cigar_mismatch,1,cigar_attr_none,3}, {cigar_match,59,cigar_attr_none,0} };
  load_vector(matches->cigar_vector,hand,3);
  load_vector(matches->cigar_vector,hand,3);
  traces[0].match_alignment.cigar_offset = 3;
  traces[0].match_alignment.cigar_length = 3;
  GT_CHECK(match_trace_get_cigar_buffer(traces+0,matches->cigar_vector)==vector_get_mem(matches->cigar_vector,cigar_element_t)+3);
  GT_CHECK_EQ_UINT(match_trace_get_cigar_length(traces+0),3);
  // Copy keeps every field that describes the match
  match_trace_t copy;
  memset(&copy,0xA5,sizeof(copy));
  traces[0].mapq_score = 37;
  traces[0].bs_strand = bs_strand_G2A;
  traces[0].strand = Reverse;
  match_trace_copy(&copy,traces+0);
  GT_CHECK(copy.sequence_name==traces[0].sequence_name);
  GT_CHECK_EQ_UINT(copy.text_position,traces[0].text_position);
  GT_CHECK_EQ_INT(copy.strand,Reverse);
  GT_CHECK_EQ_INT(copy.bs_strand,bs_strand_G2A);
  GT_CHECK_EQ_INT(copy.swg_score,traces[0].swg_score);
  GT_CHECK_EQ_UINT(copy.edit_distance,traces[0].edit_distance);
  GT_CHECK_EQ_UINT(copy.event_distance,traces[0].event_distance);
  GT_CHECK(copy.error_quality==traces[0].error_quality);
  GT_CHECK_EQ_UINT(copy.mapq_score,37);
  GT_CHECK_EQ_UINT(copy.match_alignment.cigar_offset,3);
  GT_CHECK_EQ_INT(copy.match_alignment.effective_length,traces[0].match_alignment.effective_length);
  GT_CHECK_EQ_INT(copy.type,traces[0].type);
  matches_delete(matches);
}

/*
 * Local matches wait in a pending list; moving them to the extended list
 * marks them extended
 */
GT_TEST(matches_local_pending) {
  matches_t* const matches = matches_new();
  static match_trace_t t[2];
  make_trace(t+0,50,5,5,20.0f,0,10);
  make_trace(t+1,40,6,6,20.0f,1,20);
  t[0].type = t[1].type = match_type_local;
  matches_add_match_trace_local_pending(matches,t+0);
  matches_add_match_trace_local_pending(matches,t+1);
  GT_CHECK_EQ_UINT(vector_get_used(matches->match_traces_local),2);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces_extended(matches),0);
  matches_local_pending_add_to_extended_matches(matches,NULL);
  GT_CHECK_EQ_UINT(vector_get_used(matches->match_traces_local),0);
  GT_REQUIRE(matches_get_num_match_traces_extended(matches)==2);
  match_trace_t** const ext = vector_get_mem(matches->match_traces_extended,match_trace_t*);
  GT_CHECK_EQ_INT(ext[0]->type,match_type_extended);
  GT_CHECK_EQ_INT(ext[0]->swg_score,50);
  GT_CHECK_EQ_INT(ext[1]->swg_score,40);
  GT_CHECK(ext[0]!=t+0); // a copy owned by the matches
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(matches),0);
  matches_delete(matches);
}

/*
 * Classification (hand-worked)
 */
static matches_t* matches_with(match_trace_t* const traces,const uint64_t n,const uint64_t mcs) {
  matches_t* const matches = matches_new();
  uint64_t i;
  for (i=0;i<n;++i) {
    matches_add_match_trace_insert_sorted(matches,traces+i);
    matches_counters_add(matches->counters,traces[i].edit_distance,1);
    matches_metrics_update(&matches->metrics,traces[i].event_distance,traces[i].edit_distance,traces[i].swg_score);
  }
  matches->max_complete_stratum = mcs;
  return matches;
}
GT_TEST(classify_hand) {
  static match_trace_t t[4];
  matches_t* m;
  // Unmapped
  m = matches_with(t,0,3);
  matches_classify(m);
  GT_CHECK_EQ_INT(m->classification.matches_class,matches_class_unmapped);
  GT_CHECK_EQ_INT(m->classification.delta_group,0);
  GT_CHECK_EQ_INT(m->classification.wdelta_group,0);
  matches_delete(m);
  // Unique: wdelta = mcs - best edit
  make_trace(t+0,92,2,2,30.0f,0,10);
  m = matches_with(t,1,5);
  matches_classify(m);
  GT_CHECK_EQ_INT(m->classification.matches_class,matches_class_unique);
  GT_CHECK_EQ_INT(m->classification.delta_group,-1);
  GT_CHECK_EQ_INT(m->classification.wdelta_group,3);
  matches_delete(m);
  // Multi-map with a distinct best: delta = 4-2, wdelta = min(6,4)-2
  make_trace(t+0,92,2,2,30.0f,0,10);
  make_trace(t+1,84,4,4,30.0f,0,90);
  m = matches_with(t,2,6);
  matches_classify(m);
  GT_CHECK_EQ_INT(m->classification.matches_class,matches_class_mmap);
  GT_CHECK_EQ_INT(m->classification.delta_group,2);
  GT_CHECK_EQ_INT(m->classification.wdelta_group,2);
  matches_delete(m);
  // ... mcs below the second edit distance limits wdelta
  m = matches_with(t,2,3);
  matches_classify(m);
  GT_CHECK_EQ_INT(m->classification.delta_group,2);
  GT_CHECK_EQ_INT(m->classification.wdelta_group,1);
  matches_delete(m);
  // Tie on SWG score
  make_trace(t+0,90,2,2,30.0f,0,10);
  make_trace(t+1,90,3,3,30.0f,0,90);
  m = matches_with(t,2,6);
  matches_classify(m);
  GT_CHECK_EQ_INT(m->classification.delta_group,0);
  GT_CHECK_EQ_INT(m->classification.wdelta_group,1);
  matches_delete(m);
  // Tie on edit distance with different scores
  make_trace(t+0,95,1,2,30.0f,0,10);
  make_trace(t+1,85,3,2,30.0f,0,90);
  m = matches_with(t,2,6);
  matches_classify(m);
  GT_CHECK_EQ_INT(m->classification.delta_group,0);
  matches_delete(m);
  // Primary (best score) is not the smallest edit distance: a tie
  make_trace(t+0,95,1,3,30.0f,0,10);
  make_trace(t+1,85,2,1,30.0f,0,90);
  m = matches_with(t,2,6);
  matches_classify(m);
  GT_CHECK_EQ_INT(m->classification.delta_group,0);
  GT_CHECK_EQ_INT(m->classification.wdelta_group,2); // min(6,3)-1
  matches_delete(m);
}

GT_TEST(classify_min_depth_searched) {
  static const struct { matches_class_type c; int64_t delta, wdelta; bool done; } cases[] = {
    { matches_class_unmapped, 0, 0, false },
    { matches_class_unique, -1, 2, true },
    { matches_class_unique, -1, 1, false },
    { matches_class_unique, -1, 0, false },
    { matches_class_mmap, 0, 0, true },
    { matches_class_mmap, 0, 3, true },
    { matches_class_mmap, 0, -1, false },
    { matches_class_mmap, 2, 0, false },
    { matches_class_mmap, 1, -1, false },
    { matches_class_mmap, 1, 1, true },
    { matches_class_mmap, 2, 1, false },
    { matches_class_mmap, 3, 2, true },
  };
  uint64_t i;
  for (i=0;i<sizeof(cases)/sizeof(cases[0]);++i) {
    matches_classification_t c = { cases[i].c, cases[i].delta, cases[i].wdelta };
    GT_CHECK_MSG(matches_classify_min_depth_searched(&c)==cases[i].done,"case %" PRIu64,i);
  }
}

GT_TEST(classify_max_search_error) {
  static match_trace_t t[4];
  pattern_t pattern;
  memset(&pattern,0,sizeof(pattern));
  pattern.key_length = 100;
  const uint64_t proper_length = 10; // max pattern error MAX(100/20,1) = 5
  matches_t* m;
  // Unmapped: the feasible error (mcs 0 -> 5; mcs 2 -> 3)
  m = matches_with(t,0,0);
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),5);
  m->max_complete_stratum = 2;
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),3);
  matches_delete(m);
  // Unique well separated (wdelta 3): mcs-1
  make_trace(t+0,98,1,1,30.0f,0,10);
  m = matches_with(t,1,4);
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),3);
  matches_delete(m);
  // Unique, wdelta 1: best+2, bounded by feasible 2*mcs-1
  make_trace(t+0,94,3,3,30.0f,0,10);
  m = matches_with(t,1,4);
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),5);
  m->max_complete_stratum = 3; // wdelta 0, feasible 5
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),5);
  matches_delete(m);
  // Multi-map tie, wdelta >= 0: mcs-1
  make_trace(t+0,90,2,2,30.0f,0,10);
  make_trace(t+1,90,2,2,30.0f,0,90);
  m = matches_with(t,2,4);
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),3);
  // ... wdelta < 0: best-1
  m->max_complete_stratum = 1;
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),1);
  matches_delete(m);
  // Multi-map, delta 1, wdelta <= 0: best
  make_trace(t+0,92,2,2,30.0f,0,10);
  make_trace(t+1,86,3,3,30.0f,0,90);
  m = matches_with(t,2,2);
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),2);
  // ... delta 1, wdelta 1: mcs-1
  m->max_complete_stratum = 5;
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),4);
  matches_delete(m);
  // Multi-map, delta 3: wdelta <= 0 -> best+1; wdelta 1 -> best+1; wdelta > 1 -> mcs-1
  make_trace(t+0,96,1,1,30.0f,0,10);
  make_trace(t+1,84,4,4,30.0f,0,90);
  m = matches_with(t,2,1);
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),1); // min(2, feasible 1)
  m->max_complete_stratum = 2;
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),2);
  m->max_complete_stratum = 6;
  GT_CHECK_EQ_UINT(matches_classify_compute_max_search_error(m,&pattern,proper_length),5);
  matches_delete(m);
}

GT_TEST(classify_fallbacks) {
  static match_trace_t t[200];
  search_parameters_t sp;
  search_parameters_init(&sp);
  pattern_t pattern;
  memset(&pattern,0,sizeof(pattern));
  // Viable patterns have at least one non-wildcard base
  pattern.key_length = 10; pattern.num_wildcards = 3;
  GT_CHECK(matches_classify_pattern_viable(&pattern));
  pattern.num_wildcards = 10;
  GT_CHECK(!matches_classify_pattern_viable(&pattern));
  pattern.key_length = 0; pattern.num_wildcards = 0;
  GT_CHECK(!matches_classify_pattern_viable(&pattern));
  pattern.key_length = 100; pattern.num_wildcards = 0;
  // Local alignment only for unmapped reads, unless disabled
  matches_t* m = matches_with(t,0,0);
  GT_CHECK(matches_classify_local_alignment_fallback(m,local_alignment_if_unmapped));
  GT_CHECK(!matches_classify_local_alignment_fallback(m,local_alignment_never));
  matches_delete(m);
  make_trace(t+0,100,0,0,30.0f,0,10);
  m = matches_with(t,1,0);
  GT_CHECK(!matches_classify_local_alignment_fallback(m,local_alignment_if_unmapped));
  // -s: the error is capped at best + strata after best
  GT_CHECK_EQ_UINT(matches_classify_adjust_max_error_by_strata_after_best(m,12,1),1);
  matches_delete(m);
  make_trace(t+0,92,2,2,30.0f,0,10);
  m = matches_with(t,1,0);
  GT_CHECK_EQ_UINT(matches_classify_adjust_max_error_by_strata_after_best(m,12,3),5);
  GT_CHECK_EQ_UINT(matches_classify_adjust_max_error_by_strata_after_best(m,4,3),4);
  matches_delete(m);
  m = matches_with(t,0,0);
  GT_CHECK_EQ_UINT(matches_classify_adjust_max_error_by_strata_after_best(m,7,1),7);
  matches_delete(m);
  // Max matches searched: once max_searched_matches matches are found, the
  // search stops when the edit-distance bound of the last one (score 93 of
  // 100 with penalties 1,4,6,1: at most 2 mismatches or one 1-base gap, so
  // 2) is within the max complete stratum
  sp.select_parameters.max_searched_matches = 5;
  sp.select_parameters.min_reported_strata_nominal = 0;
  uint64_t i;
  for (i=0;i<5;++i) make_trace(t+i,93,1,1,30.0f,0,10+i);
  m = matches_with(t,4,2);
  GT_CHECK(!matches_classify_max_matches_searched(m,&sp,&pattern)); // only 4
  matches_delete(m);
  m = matches_with(t,5,2);
  GT_CHECK(matches_classify_max_matches_searched(m,&sp,&pattern));
  // These ties are deep enough already: no neighbourhood fallback
  GT_CHECK(!matches_classify_neighbourhood_fallback(m,&sp,&pattern));
  sp.select_parameters.min_reported_strata_nominal = 1; // -m disables the limit
  GT_CHECK(!matches_classify_max_matches_searched(m,&sp,&pattern));
  sp.select_parameters.min_reported_strata_nominal = 0;
  m->max_complete_stratum = 1;
  GT_CHECK(!matches_classify_max_matches_searched(m,&sp,&pattern));
  matches_delete(m);
  // A poor fifth match leaves room for better ones unless mcs covers its bound
  for (i=0;i<5;++i) make_trace(t+i,40,10,10,30.0f,0,10+i);
  m = matches_with(t,5,1);
  GT_CHECK(!matches_classify_max_matches_searched(m,&sp,&pattern));
  m->max_complete_stratum = 100;
  GT_CHECK(matches_classify_max_matches_searched(m,&sp,&pattern));
  matches_delete(m);
  // Unmapped and shallow: fall back to neighbourhood search
  m = matches_with(t,0,0);
  GT_CHECK(matches_classify_neighbourhood_fallback(m,&sp,&pattern));
  matches_delete(m);
  // Unique and well separated: done
  make_trace(t+0,100,0,0,30.0f,0,10);
  m = matches_with(t,1,3);
  GT_CHECK(!matches_classify_neighbourhood_fallback(m,&sp,&pattern));
  matches_delete(m);
}

/*
 * Predictors normalise distances by the read length and scores by the
 * maximum score (match score x read length)
 */
GT_TEST(predictors) {
  matches_metrics_t metrics;
  matches_metrics_init(&metrics);
  matches_metrics_set_read_length(&metrics,100);
  matches_metrics_set_swg_match_score(&metrics,2);
  matches_metrics_set_proper_length(&metrics,8.0);
  matches_metrics_set_region_profile_metrics(&metrics,24,40,0.5);
  matches_metrics_update(&metrics,3,4,170);
  matches_metrics_update(&metrics,6,8,120);
  matches_predictors_t p;
  memset(&p,0xff,sizeof(p));
  matches_predictors_compute_unique(&p,&metrics,4,3,170,31.5,5.0);
  GT_CHECK(fabs(p.primary_edit_distance_norm-0.96) < 1e-12);
  GT_CHECK(fabs(p.primary_event_distance_norm-0.97) < 1e-12);
  GT_CHECK(fabs(p.primary_swg_score_norm-0.85) < 1e-12);
  GT_CHECK(fabs(p.primary_error_quality-31.5) < 1e-12);
  GT_CHECK(fabs(p.primary_template_size_sigma_norm-0.5) < 1e-12); // 5 of MAX_TEMPLATE_LENGTH_SIGMAS
  GT_CHECK(p.subdominant_edit_distance_norm==0.0 && p.subdominant_event_distance_norm==0.0);
  GT_CHECK(p.subdominant_swg_score_norm==0.0 && p.subdominant_template_size_sigma_norm==1.0);
  matches_predictors_compute_mmaps(&p,&metrics,4,3,170,31.5,0.0);
  GT_CHECK(fabs(p.subdominant_edit_distance_norm-0.92) < 1e-12);
  GT_CHECK(fabs(p.subdominant_event_distance_norm-0.94) < 1e-12);
  GT_CHECK(fabs(p.subdominant_swg_score_norm-0.60) < 1e-12);
  GT_CHECK(p.subdominant_template_size_sigma_norm==0.0);
  GT_CHECK(p.primary_template_size_sigma_norm==0.0);
  // Single-end: mappability and search scope come from the matches
  static match_trace_t t[3];
  make_trace(t+0,170,3,4,31.5f,0,10);
  make_trace(t+1,120,6,8,25.0f,0,90);
  make_trace(t+2,110,7,8,25.0f,0,190);
  matches_t* const m = matches_with(t,3,9);
  m->metrics.read_length = 100;
  m->metrics.swg_match_score = 2;
  m->metrics.proper_length = 8.0;
  m->metrics.avg_region_length = 24;
  m->metrics.max_region_length = 40;
  m->metrics.kmer_frequency = 0.5;
  m->metrics.candidates_accepted = 11;
  matches_predictors_compute_se(&p,m);
  GT_CHECK_EQ_INT(m->classification.matches_class,matches_class_mmap);
  GT_CHECK_EQ_UINT(p.candidates_subdominant,2); // two at edit distance 8
  GT_CHECK_EQ_UINT(p.candidates_accepted,11);
  GT_CHECK_EQ_UINT(p.mcs_end1,9);
  GT_CHECK_EQ_UINT(p.mcs_end2,0);
  GT_CHECK(fabs(p.primary_template_size_sigma_norm-1.0) < 1e-12);
  GT_CHECK(fabs(p.sequence_length_norm-12.5) < 1e-12);
  GT_CHECK(fabs(p.avg_region_length_norm-3.0) < 1e-12);
  GT_CHECK(fabs(p.max_region_length_norm-5.0) < 1e-12);
  GT_CHECK(fabs(p.kmer_frequency-0.5) < 1e-12);
  GT_CHECK(fabs(p.subdominant_edit_distance_norm-0.92) < 1e-12);
  GT_CHECK_EQ_UINT(p.mapq_end1,0);
  GT_CHECK_EQ_UINT(p.mapq_end2,0);
  matches_delete(m);
}

/*
 * Selection: -M keeps the first M matches; -m keeps whole strata
 */
GT_TEST(select_parameters) {
  select_parameters_t s;
  memset(&s,0xff,sizeof(s));
  select_parameters_init(&s);
  GT_CHECK(s.min_reported_strata==0.0);
  GT_CHECK_EQ_UINT(s.max_reported_matches,5);
  GT_CHECK_EQ_UINT(s.max_searched_matches,100);
  GT_CHECK_EQ_UINT(s.max_searched_paired_matches,20);
  // -m: a fraction is of the read length, 1 or more is a count
  s.min_reported_strata = 0.05;
  select_parameters_instantiate_values(&s,200);
  GT_CHECK_EQ_UINT(s.min_reported_strata_nominal,10);
  s.min_reported_strata = 3.0;
  select_parameters_instantiate_values(&s,200);
  GT_CHECK_EQ_UINT(s.min_reported_strata_nominal,3);
  s.min_reported_strata = 0.0;
  select_parameters_instantiate_values(&s,200);
  GT_CHECK_EQ_UINT(s.min_reported_strata_nominal,0);
  // SE searches look for as many matches as -M reports; PE for fixed limits
  archive_search_t* const as = calloc(1,sizeof(archive_search_t));
  select_parameters_init(&as->search_parameters.select_parameters);
  as->search_parameters.select_parameters.max_reported_matches = 7;
  archive_select_configure_se(as);
  GT_CHECK_EQ_UINT(as->search_parameters.select_parameters.max_searched_matches,7);
  GT_CHECK_EQ_UINT(as->search_parameters.select_parameters.max_searched_paired_matches,0);
  archive_select_configure_pe(as);
  GT_CHECK_EQ_UINT(as->search_parameters.select_parameters.max_searched_matches,100);
  GT_CHECK_EQ_UINT(as->search_parameters.select_parameters.max_searched_paired_matches,20);
  GT_CHECK_EQ_UINT(as->search_parameters.select_parameters.max_reported_matches,7);
  free(as);
}

GT_TEST(select_se_matches) {
  static match_trace_t t[16];
  // Edit distances 0, 1, 1, 2, 2, 2, 5 (sorted by score)
  static const uint64_t edits[] = {0,1,1,2,2,2,5};
  uint64_t i;
  for (i=0;i<7;++i) make_trace(t+i,100-(int32_t)(10*edits[i])-(int32_t)i,edits[i],edits[i],30.0f,0,10+i);
  select_parameters_t s;
  select_parameters_init(&s);
  // -M 4: the best four
  matches_t* m = matches_with(t,7,0);
  s.max_reported_matches = 4;
  s.min_reported_strata_nominal = 0;
  archive_select_se_matches(&s,m);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(m),4);
  GT_CHECK(matches_get_match_traces(m)[3]==t+3);
  matches_delete(m);
  // -M larger than the number of matches: unchanged
  m = matches_with(t,7,0);
  s.max_reported_matches = 100;
  archive_select_se_matches(&s,m);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(m),7);
  matches_delete(m);
  // -m 1 (with -M 2): the whole first stratum and the next ones within -M;
  // matches up to two strata beyond those are kept (base behaviour)
  m = matches_with(t,7,0);
  s.max_reported_matches = 2;
  s.min_reported_strata_nominal = 1;
  archive_select_se_matches(&s,m);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(m),6); // 0:1 taken; 1:2 would exceed 2 -> strata 1, keep edit <= 2
  for (i=0;i<matches_get_num_match_traces(m);++i) GT_CHECK(matches_get_match_traces(m)[i]->edit_distance<=2);
  matches_delete(m);
  // Discard keeps edit distance <= strata+1, in order
  m = matches_with(t,7,0);
  archive_select_se_matches_discard(m,0);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(m),3);
  GT_CHECK(matches_get_match_traces(m)[0]==t+0 && matches_get_match_traces(m)[2]==t+2);
  archive_select_se_matches_discard(m,10);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(m),3);
  matches_delete(m);
  // -m on unmapped reads leaves them unmapped
  m = matches_with(t,0,0);
  archive_select_se_matches(&s,m);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(m),0);
  matches_delete(m);
}

/*
 * MAPQ: probability to Phred, and the logistic model
 */
GT_TEST(score_probability_and_logit) {
  GT_CHECK_EQ_UINT(archive_score_probability_to_mapq(0.0),0);
  GT_CHECK_EQ_UINT(archive_score_probability_to_mapq(0.5),3);      // 3.01
  GT_CHECK_EQ_UINT(archive_score_probability_to_mapq(0.9),10);
  GT_CHECK_EQ_UINT(archive_score_probability_to_mapq(0.99),20);
  GT_CHECK_EQ_UINT(archive_score_probability_to_mapq(0.999),30);
  GT_CHECK_EQ_UINT(archive_score_probability_to_mapq(0.9999),40);
  GT_CHECK_EQ_UINT(archive_score_probability_to_mapq(0.75),6);     // 6.02
  gt_rng_t rng = gt_rng_new(112);
  uint64_t trial;
  for (trial=0;trial<200;++trial) {
    archive_score_logit_coeff_t c;
    matches_predictors_t p;
    matches_classification_t cl;
    double* const cc = (double*)&c;
    uint64_t k;
    for (k=0;k<sizeof(c)/sizeof(double);++k) cc[k] = ((double)gt_rng_below(&rng,2001)-1000.0)/1000.0;
    p.primary_edit_distance_norm = gt_rng_below(&rng,100)/100.0;
    p.primary_event_distance_norm = gt_rng_below(&rng,100)/100.0;
    p.primary_swg_score_norm = gt_rng_below(&rng,100)/100.0;
    p.primary_error_quality = gt_rng_below(&rng,40);
    p.primary_template_size_sigma_norm = gt_rng_below(&rng,100)/100.0;
    p.mapq_end1 = (uint8_t)gt_rng_below(&rng,61);
    p.mapq_end2 = (uint8_t)gt_rng_below(&rng,61);
    p.subdominant_edit_distance_norm = gt_rng_below(&rng,100)/100.0;
    p.subdominant_event_distance_norm = gt_rng_below(&rng,100)/100.0;
    p.subdominant_swg_score_norm = gt_rng_below(&rng,100)/100.0;
    p.subdominant_template_size_sigma_norm = gt_rng_below(&rng,100)/100.0;
    p.candidates_subdominant = gt_rng_below(&rng,10);
    p.candidates_accepted = gt_rng_below(&rng,10);
    p.mcs_end1 = gt_rng_below(&rng,10);
    p.mcs_end2 = gt_rng_below(&rng,10);
    p.sequence_length_norm = gt_rng_below(&rng,100)/10.0;
    p.avg_region_length_norm = gt_rng_below(&rng,100)/10.0;
    p.max_region_length_norm = gt_rng_below(&rng,100)/10.0;
    p.kmer_frequency = gt_rng_below(&rng,100)/100.0;
    cl.matches_class = matches_class_mmap;
    cl.delta_group = (int64_t)gt_rng_below(&rng,5)-1;
    cl.wdelta_group = (int64_t)gt_rng_below(&rng,5)-1;
    const double x = c.coeff_intercept
        + c.coeff_primary_edit_distance_norm*p.primary_edit_distance_norm
        + c.coeff_primary_event_distance_norm*p.primary_event_distance_norm
        + c.coeff_primary_swg_score_norm*p.primary_swg_score_norm
        + c.coeff_primary_error_quality*p.primary_error_quality
        + c.coeff_primary_template_size_sigma_norm*p.primary_template_size_sigma_norm
        + c.coeff_mapq_end1*p.mapq_end1 + c.coeff_mapq_end2*p.mapq_end2
        + c.coeff_subdominant_edit_distance_norm*p.subdominant_edit_distance_norm
        + c.coeff_subdominant_event_distance_norm*p.subdominant_event_distance_norm
        + c.coeff_subdominant_swg_score_norm*p.subdominant_swg_score_norm
        + c.coeff_subdominant_template_size_sigma_norm*p.subdominant_template_size_sigma_norm
        + c.coeff_candidates_subdominant*(double)p.candidates_subdominant
        + c.coeff_candidates_accepted*(double)p.candidates_accepted
        + c.coeff_mcs_end1*(double)p.mcs_end1 + c.coeff_mcs_end2*(double)p.mcs_end2
        + c.coeff_delta_group*(double)cl.delta_group + c.coeff_wdelta_group*(double)cl.wdelta_group
        + c.coeff_sequence_length_norm*p.sequence_length_norm
        + c.coeff_avg_region_length_norm*p.avg_region_length_norm
        + c.coeff_max_region_length_norm*p.max_region_length_norm
        + c.coeff_kmer_frequency*p.kmer_frequency;
    const double expected = 1.0/(1.0+exp(-x));
    const double got = archive_score_logit(&c,&p,&cl);
    GT_CHECK_MSG(fabs(got-expected) < 1e-9,"trial %" PRIu64 ": %.12f != %.12f",trial,got,expected);
    // unique/mmap pick their coefficient set
    archive_score_logit_model_t model;
    model.unique_logit_coeff = c;
    memset(&model.mmap_logit_coeff,0,sizeof(model.mmap_logit_coeff));
    GT_CHECK(fabs(archive_score_logit_unique(&model,&p,&cl)-expected) < 1e-9);
    GT_CHECK(fabs(archive_score_logit_mmap(&model,&p,&cl)-0.5) < 1e-12);
  }
}

/*
 * Single-end MAPQ models on hand-built matches
 */
static archive_search_t* new_archive_search(const mapq_model_t model,const mapping_mode_t mode) {
  archive_search_t* const as = calloc(1,sizeof(archive_search_t));
  search_parameters_init(&as->search_parameters);
  as->search_parameters.mapq_model_se = model;
  as->search_parameters.mapping_mode = mode;
  return as;
}
static void set_mappability(matches_t* const m) {
  m->metrics.read_length = 100;
  m->metrics.swg_match_score = 1;
  m->metrics.proper_length = 9.0;
  m->metrics.avg_region_length = 20;
  m->metrics.max_region_length = 30;
  m->metrics.kmer_frequency = 0.5;
}
static uint8_t score_se(
    const mapq_model_t model,const mapping_mode_t mode,match_trace_t* const t,const uint64_t n,
    const uint64_t mcs,const match_type type) {
  archive_search_t* const as = new_archive_search(model,mode);
  matches_t* const m = matches_with(t,n,mcs);
  set_mappability(m);
  matches_get_primary_match(m)->type = type;
  matches_get_primary_match(m)->mapq_score = 77; // marker: unchanged by 'none'
  archive_score_matches_se(as,m);
  const uint8_t mapq = matches_get_primary_match(m)->mapq_score;
  matches_get_primary_match(m)->type = match_type_regular;
  matches_delete(m);
  free(as);
  return mapq;
}
// Logistic MAPQ expected from the model's published coefficients
static uint8_t expected_logit_mapq(const archive_score_logit_coeff_t* const c,match_trace_t* const t,const uint64_t n,const uint64_t mcs) {
  matches_t* const m = matches_with(t,n,mcs);
  set_mappability(m);
  matches_classify(m);
  matches_predictors_t p;
  matches_predictors_compute_se(&p,m);
  const double pr = archive_score_logit(c,&p,&m->classification);
  matches_delete(m);
  const double q = round(-10.0*log10(1.0-pr));
  return (uint8_t)(q < 59.0 ? q : 59.0);
}
GT_TEST(score_se_models) {
  static match_trace_t t[4];
  // Exact unique read well separated (wdelta 3): 60 (gem), 243 (classify)
  make_trace(t+0,100,0,0,30.0f,0,10);
  GT_CHECK_EQ_UINT(score_se(mapq_model_gem,mapping_adaptive_filtering_fast,t,1,3,match_type_regular),60);
  GT_CHECK_EQ_UINT(score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,1,3,match_type_regular),243);
  GT_CHECK_EQ_UINT(score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,1,12,match_type_regular),250);
  GT_CHECK_EQ_UINT(score_se(mapq_model_none,mapping_adaptive_filtering_fast,t,1,3,match_type_regular),77);
  // Unique, wdelta <= 1: the logistic model of the mapping mode, capped at 59
  make_trace(t+0,88,3,3,30.0f,0,10);
  const uint8_t fast = score_se(mapq_model_gem,mapping_adaptive_filtering_fast,t,1,4,match_type_regular);
  const uint8_t sens = score_se(mapq_model_gem,mapping_hybrid_sensitive,t,1,4,match_type_regular);
  GT_CHECK_EQ_UINT(fast,expected_logit_mapq(&logit_model_single_end_fast.unique_logit_coeff,t,1,4));
  GT_CHECK_EQ_UINT(sens,expected_logit_mapq(&logit_model_single_end_sensitive.unique_logit_coeff,t,1,4));
  GT_CHECK(fast<=59 && sens<=59);
  {
    search_parameters_t sp;
    search_parameters_init(&sp);
    sp.mapping_mode = mapping_adaptive_filtering_fast;
    GT_CHECK(archive_score_matches_se_get_model(&sp)==&logit_model_single_end_fast);
    sp.mapping_mode = mapping_hybrid_sensitive;
    GT_CHECK(archive_score_matches_se_get_model(&sp)==&logit_model_single_end_sensitive);
    sp.mapping_mode = mapping_neighborhood_search_brute_force;
    GT_CHECK(archive_score_matches_se_get_model(&sp)==&logit_model_single_end_sensitive);
  }
  // Classify: unique with wdelta <= 1 in [110,170]
  const uint8_t cu = score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,1,4,match_type_regular);
  GT_CHECK(cu>=110 && cu<=170);
  // Ties: 0 (gem); classify 0 or 2
  make_trace(t+0,100,0,0,30.0f,0,10);
  make_trace(t+1,100,0,0,30.0f,0,90);
  GT_CHECK_EQ_UINT(score_se(mapq_model_gem,mapping_adaptive_filtering_fast,t,2,3,match_type_regular),0);
  GT_CHECK_EQ_UINT(score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,2,3,match_type_regular),0);
  // Classify tie of equal events but different edit and score: 2
  make_trace(t+0,95,1,1,30.0f,0,10);
  make_trace(t+1,90,1,2,30.0f,0,90);
  // (delta 0 because the events tie)
  GT_CHECK_EQ_UINT(score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,2,6,match_type_regular),2);
  // Multi-map with a clear best (delta 3, wdelta 3): 60 (gem), 223 (classify)
  make_trace(t+0,100,0,0,30.0f,0,10);
  make_trace(t+1,85,3,3,30.0f,0,90);
  GT_CHECK_EQ_UINT(score_se(mapq_model_gem,mapping_adaptive_filtering_fast,t,2,5,match_type_regular),60);
  GT_CHECK_EQ_UINT(score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,2,5,match_type_regular),223);
  make_trace(t+1,40,12,12,30.0f,0,90);
  GT_CHECK_EQ_UINT(score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,2,20,match_type_regular),230);
  // ... wdelta 1: logistic mmap model (gem, <= 59), classify in [10,70]
  make_trace(t+1,85,3,3,30.0f,0,90);
  const uint8_t mm = score_se(mapq_model_gem,mapping_adaptive_filtering_fast,t,2,1,match_type_regular);
  GT_CHECK_EQ_UINT(mm,expected_logit_mapq(&logit_model_single_end_fast.mmap_logit_coeff,t,2,1));
  const uint8_t cm = score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,2,1,match_type_regular);
  GT_CHECK(cm>=10 && cm<=70);
  // Local primary: 1 when short (gem: effective length < 40) or always (classify)
  make_trace(t+0,30,0,0,30.0f,0,10);
  t[0].match_alignment.effective_length = 39;
  GT_CHECK_EQ_UINT(score_se(mapq_model_gem,mapping_adaptive_filtering_fast,t,1,3,match_type_local),1);
  t[0].match_alignment.effective_length = 40;
  GT_CHECK_EQ_UINT(score_se(mapq_model_gem,mapping_adaptive_filtering_fast,t,1,3,match_type_local),60);
  GT_CHECK_EQ_UINT(score_se(mapq_model_classify,mapping_adaptive_filtering_fast,t,1,3,match_type_local),1);
  // No alignment model: no MAPQ
  {
    archive_search_t* const as = new_archive_search(mapq_model_gem,mapping_adaptive_filtering_fast);
    as->search_parameters.match_alignment_model = match_alignment_model_none;
    make_trace(t+0,100,0,0,30.0f,0,10);
    matches_t* const m = matches_with(t,1,3);
    set_mappability(m);
    matches_get_primary_match(m)->mapq_score = 77;
    archive_score_matches_se(as,m);
    GT_CHECK_EQ_UINT(matches_get_primary_match(m)->mapq_score,77);
    matches_delete(m);
    free(as);
  }
  // The logistic helper caps at 59: perfect alignment, mcs 4 gives a linear
  // predictor near 18 in the fast unique model (Phred 78 before the cap)
  {
    search_parameters_t sp;
    search_parameters_init(&sp);
    matches_predictors_t p;
    memset(&p,0,sizeof(p));
    p.primary_swg_score_norm = 1.0;
    p.primary_event_distance_norm = 1.0;
    p.primary_edit_distance_norm = 1.0;
    p.mcs_end1 = 4;
    matches_classification_t cl = { matches_class_unique, -1, 1 };
    GT_CHECK_EQ_UINT(archive_score_matches_se_logit_unique(&sp,&p,&cl),59);
    p.mcs_end1 = 1; // linear predictor near -1: Phred 1
    GT_CHECK_EQ_UINT(archive_score_matches_se_logit_unique(&sp,&p,&cl),
        (uint8_t)round(-10.0*log10(1.0-archive_score_logit_unique(&logit_model_single_end_fast,&p,&cl))));
  }
}

GT_CASES(
  GT_CASE(cigar_append_merges_runs),
  GT_CASE(cigar_append_starts_new_cigar),
  GT_CASE(cigar_buffer_sentinel),
  GT_CASE(cigar_insert_element),
  GT_CASE(cigar_distances),
  GT_CASE(cigar_error_quality),
  GT_CASE(cigar_print),
  GT_CASE(counters_random_ops),
  GT_CASE(counters_matches_to_report),
  GT_CASE(metrics_update),
  GT_CASE(matches_insert_sorted),
  GT_CASE(matches_preserve_sorted),
  GT_CASE(matches_sort_by_genomic_position),
  GT_CASE(matches_index_lookup),
  GT_CASE(matches_strata_bookkeeping),
  GT_CASE(matches_recompute_and_accessors),
  GT_CASE(matches_local_pending),
  GT_CASE(classify_hand),
  GT_CASE(classify_min_depth_searched),
  GT_CASE(classify_max_search_error),
  GT_CASE(classify_fallbacks),
  GT_CASE(predictors),
  GT_CASE(select_parameters),
  GT_CASE(select_se_matches),
  GT_CASE(score_probability_and_logit),
  GT_CASE(score_se_models)
)
