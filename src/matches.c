/*
 * PROJECT: GEMMapper
 * FILE: matches.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure to store alignment matches {sequence,position,strand,CIGAR}
 */

#include "matches.h"

/*
 * Constants
 */
#define MATCHES_INIT_COUNTERS          200
#define MATCHES_INIT_INTERVAL_MATCHES  100
#define MATCHES_INIT_GLOBAL_MATCHES   1000
#define MATCHES_INIT_CIGAR_OPS        5000

/*
 * Setup
 */
GEM_INLINE matches_t* matches_new() {
  // Allocate handler
  matches_t* const matches = mm_alloc(matches_t);
  // Text Collection Buffer
  matches->text_collection = text_collection_new();
  // Matches Counters
  matches->counters = vector_new(MATCHES_INIT_COUNTERS,uint64_t);
  // Interval Matches
  matches->interval_matches = vector_new(MATCHES_INIT_INTERVAL_MATCHES,match_interval_t);
  // Position Matches
  matches->global_matches = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t);
  matches->start_gmatches = ihash_new();
  // CIGAR buffer
  matches->cigar_buffer = vector_new(MATCHES_INIT_CIGAR_OPS,cigar_element_t);
  // Restore Point (RP)
  matches->rp_counters = vector_new(MATCHES_INIT_COUNTERS,uint64_t);
  matches->rp_interval_matches_used = 0;
  matches->rp_global_matches_used = 0;
  matches->rp_cigar_buffer_used = 0;
  // Pre-computed Data
  matches->total_matches = 0;
  matches->last_computed_interval_matches_used = 0;
  matches->last_computed_global_matches_used = 0;
  // Return
  return matches;
}
GEM_INLINE void matches_clear(matches_t* const matches) {
  MATCHES_CHECK(matches);
  text_collection_clear(matches->text_collection);
  vector_clear(matches->counters);
  vector_clear(matches->interval_matches);
  vector_clear(matches->global_matches);
  ihash_clear(matches->start_gmatches);
  vector_clear(matches->cigar_buffer);
  matches->total_matches = 0;
  matches->last_computed_interval_matches_used = 0;
  matches->last_computed_global_matches_used = 0;
}
GEM_INLINE void matches_delete(matches_t* const matches) {
  MATCHES_CHECK(matches);
  // Delete all
  text_collection_delete(matches->text_collection);
  vector_delete(matches->counters);
  vector_delete(matches->interval_matches);
  vector_delete(matches->global_matches);
  ihash_delete(matches->start_gmatches);
  vector_delete(matches->cigar_buffer);
  // Delete handler
  mm_free(matches);
}
/*
 * Counters
 */
GEM_INLINE uint64_t matches_get_num_matches(matches_t* const matches) {
  // Check if pre-computed value still holds
  const uint64_t interval_matches_used = vector_get_used(matches->interval_matches);
  const uint64_t global_matches_used = vector_get_used(matches->global_matches);
  if (matches->last_computed_interval_matches_used!=interval_matches_used ||
      matches->last_computed_global_matches_used!=global_matches_used) {
    // Calculate again the number of matches
    matches->last_computed_interval_matches_used = 0;
    matches->last_computed_global_matches_used = 0;
    const uint64_t num_counters = vector_get_used(matches->counters);
    const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
    uint64_t i = 0, acc = 0;
    while (i < num_counters) {
      acc += counters[i++];
    }
    matches->total_matches = acc;
  }
  return matches->total_matches;
}
GEM_INLINE uint64_t matches_counters_compact(matches_t* const matches) {
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  int64_t i = vector_get_used(matches->counters)-1;
  while (i>=0 && counters[i]==0) --i;
  vector_set_used(matches->counters,++i);
  return i;
}
GEM_INLINE uint64_t matches_counters_get_min_matching_stratum(matches_t* const matches) {
  const uint64_t num_counters = vector_get_used(matches->counters);
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  uint64_t i = 0;
  while (i<num_counters && counters[i]==0) ++i;
  return i+1;
}
/*
 * Trace-Matches
 */
GEM_INLINE cigar_element_t* match_trace_get_cigar_array(const matches_t* const matches,const match_trace_t* const match_trace) {
  return vector_get_elm(matches->cigar_buffer,match_trace->cigar_buffer_offset,cigar_element_t);
}
GEM_INLINE uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace) {
  return match_trace->cigar_length;
}
/*
 * Adding Matches
 */
GEM_INLINE void matches_counters_add(
    matches_t* const matches,const uint64_t distance,const uint64_t num_matches) {
  vector_t* const counters = matches->counters;
  // Reserve Memory
  if (distance >= vector_get_used(counters)) {
    vector_reserve(counters,distance+1,true);
    vector_set_used(counters,distance+1);
  }
  // Add matches
  *vector_get_elm(counters,distance,uint64_t) += num_matches;
}
GEM_INLINE void matches_index_match(matches_t* const matches,match_trace_t* const match_trace) {
  // Store begin position of the match (as to fast index matches)
  ihash_insert(matches->start_gmatches,match_trace->position,match_trace);
}
GEM_INLINE match_trace_t* matches_lookup_match(matches_t* const matches,const uint64_t position) {
  return ihash_get(matches->start_gmatches,position,match_trace_t);
}
GEM_INLINE void matches_add_match_trace_t(
    matches_t* const matches,match_trace_t* const match_trace) {
  vector_insert(matches->global_matches,*match_trace,match_trace_t);
}
GEM_INLINE void matches_add_match_trace_(
    matches_t* const matches,const uint64_t trace_offset,
    const uint64_t position,const uint64_t distance,const strand_t strand) {
  match_trace_t* match_trace;
  vector_alloc_new(matches->global_matches,match_trace_t,match_trace);
  match_trace->position = position;
  match_trace->trace_offset = trace_offset;
  match_trace->distance = distance;
  match_trace->strand = strand;
  match_trace->score = 0;
}
GEM_INLINE void matches_add_interval_match(
    matches_t* const matches,
    const uint64_t hi,const uint64_t lo,
    const uint64_t length,const uint64_t distance,const strand_t strand) {
  // Check non-empty interval
  const uint64_t num_matches = hi-lo;
  if (gem_expect_false(num_matches==0)) return;
  // Setup the interval-match
  match_interval_t match_interval = { hi, lo, NULL, length, distance, strand };
  // Update counters
  matches_counters_add(matches,distance,num_matches);
  // Add the interval
  vector_insert(matches->interval_matches,match_interval,match_interval_t);
}
GEM_INLINE void matches_add_interval_set(
    matches_t* const matches,interval_set_t* const interval_set) {
  // TODO NOT_IMPLEMENTED
  GEM_NOT_IMPLEMENTED();
  // TODO
}
/*
 * Handling matches
 */
GEM_INLINE void matches_reverse_CIGAR(
    matches_t* const matches,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  const uint64_t middle_point = cigar_length/2;
  uint64_t i;
  for (i=0;i<middle_point;++i) {
    cigar_element_t* const origin = cigar_buffer + i;
    cigar_element_t* const flipped = cigar_buffer + (cigar_length-1-i);
    SWAP(*origin,*flipped);
    if (origin->type == cigar_mismatch) origin->length = dna_complement(origin->length);
    if (origin->type == cigar_mismatch) origin->length = dna_complement(origin->length);
  }
  if (cigar_length%2) {
    cigar_element_t* const middle = cigar_buffer + middle_point;
    if (middle->type == cigar_mismatch) middle->length = dna_complement(middle->length);
  }
}
GEM_INLINE void matches_reverse_CIGAR_colorspace(
    matches_t* const matches,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  const uint64_t middle_point = cigar_length/2;
  uint64_t i;
  for (i=0;i<middle_point;++i) {
    cigar_element_t* const origin = cigar_buffer + i;
    cigar_element_t* const flipped = cigar_buffer + (cigar_length-1-i);
    SWAP(*origin,*flipped);
  }
}
GEM_INLINE uint64_t matches_get_effective_length(
    matches_t* const matches,const uint64_t read_length,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  const cigar_element_t* cigar_element = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  int64_t i, effective_length = read_length;
  for (i=0;i<cigar_length;++i,++cigar_element) {
    switch (cigar_element->type) {
      case cigar_ins:
        effective_length += cigar_element->length;
        break;
      case cigar_del:
      case cigar_soft_trim:
        effective_length -= cigar_element->length;
        break;
      default:
        break;
    }
  }
  GEM_INTERNAL_CHECK(effective_length<0,"Match effective length is negative");
  return effective_length;
}
/*
 * Sorting Matches
 */
int match_trace_cmp_distance(const match_trace_t* const a,const match_trace_t* const b) {
  return a->distance - b->distance;
}
GEM_INLINE void matches_sort_by_distance(matches_t* const matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(matches->global_matches,match_trace_t),
      vector_get_used(matches->global_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_distance);
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//
//#define FMI_MATCHES_INIT_AUG_MISMATCHES 3
//#define FMI_MATCHES_INIT_INTERVALS 10000
//#define FMI_MATCHES_INIT_POS 10000
//#define FMI_MATCHES_INIT_DEC 100
//#define FMI_MATCHES_NUM_QUERIES 100
//#define FMI_MULTIMATCHES_INIT_GROUPS 2
//#define FMI_MULTIMATCHES_INIT_MATCHES 500
//// Sorting instance (sort filter queries)
//#define GENSORT_NAME filter_query_t_sort
//#define GENSORT_INTTYPE idx_t
//#define GENSORT_TYPE filter_query
//#define GENSORT_KEYTYPE idx_t
//#define GENSORT_GETKEY(a) a.end_key_pos
//#define GENSORT_COMPAREKEYS(k1,k2) (k1<k2)
//#define GENSORT_USEPOINTERS
//#include "gensort.h"
//// Sorting instance (sort idx positions)
//#define GENSORT_NAME idx_t_sort
//#define GENSORT_INTTYPE idx_t
//#define GENSORT_TYPE idx_t
//#define GENSORT_KEYTYPE idx_t
//#define GENSORT_GETKEY(a) a
//#define GENSORT_COMPAREKEYS(k1,k2) (k1<k2)
//#define GENSORT_USEPOINTERS
//#include "gensort.h"
//
////
//// Private functions
////
//GEM_INLINE
//void matches_prepare(matches* const matches,const uint64_t max_mismatches) {
//  // We make room for the counters
//  register const uint64_t top_mismatches=max_mismatches+1;
//  register const uint64_t old_top_mismatches=vector_get_used(matches->rbuf_counts);
//  register uint64_t i;
//  if (top_mismatches>old_top_mismatches) {
//    vector_reserve(matches->rbuf_counts,top_mismatches);
//    register uint64_t* cntrs=vector_get_mem(matches->rbuf_counts);
//    for (i=old_top_mismatches;i<top_mismatches;++i) cntrs[i]=0ull;
//    vector_set_used(matches->rbuf_counts,top_mismatches);
//  }
//  return;
//}
//GEM_INLINE void fm_matches_append_search_results(
//    fmi_search_parameters* const search_params,
//    matches* matches, vector* found, const uint64_t num_wildcards) {
//  matches_prepare(matches, search_params->max_mismatches);
//  register uint64_t* cntrs=vector_get_mem(matches->rbuf_counts);
//
//  // Update the counters value
//  INTERVAL_ITERATE(found) {
//    register const uint64_t num_matches_interval = interval->hi-interval->lo;
//    cntrs[interval->misms+num_wildcards] += num_matches_interval;
//    matches->num_int_matches += num_matches_interval;
//  } END_INTERVAL_ITERATE;
//
//  // We append the intervals which are present in the buffer
//  vector_reserve(matches->rbuf_int,vector_get_used(matches->rbuf_int) + vector_get_used(found));
//  register int_match_t* out = vector_get_mem(matches->rbuf_int);
//  out+=vector_get_used(matches->rbuf_int);
//  INTERVAL_ITERATE(found) {
//    out->mismatches=interval->misms+num_wildcards;
//    out->lo=interval->lo;
//    out->hi=interval->hi;
//    out->key_id=search_params->key_id;
//    ++out;
//  } END_INTERVAL_ITERATE;
//  vector_set_used(matches->rbuf_int, out-((int_match_t*)vector_get_mem(matches->rbuf_int)));
//}
//GEM_INLINE void fm_matches_update_max_complete_stratum(
//    matches* matches,const int64_t max_complete_stratum) {
//  matches->max_complete_stratum = UNSAFE_MIN(matches->max_complete_stratum,max_complete_stratum);
//}
//GEM_INLINE uint64_t fm_matches_add_search_key(
//    const ch_t* key,const uint64_t key_len,const direction_t search_direction,matches* const matches) {
//  register keys_info* key_info;
//  register ch_t* key_buffer;
//  // Assign next key
//  register const uint64_t key_used=vector_get_used(matches->qbuf_keys_info);
//  // Allocate memory and store the key
//  vector_reserve(matches->qbuf_keys_info,key_used+1);
//  key_info=vector_get_mem(matches->qbuf_keys_info);
//  key_info+=key_used;
//  register const uint64_t char_used=vector_get_used(matches->qbuf_keys_buffer);
//  key_info->displacement=char_used;
//  key_info->direction=search_direction;
//  key_info->key_len=key_len;
//  vector_inc_used(matches->qbuf_keys_info);
//  vector_reserve(matches->qbuf_keys_buffer,char_used+key_len+1);
//  key_buffer=vector_get_mem(matches->qbuf_keys_buffer);
//  key_buffer+=char_used;
//  memcpy(key_buffer,key,key_len*sizeof(ch_t));
//  key_buffer[key_len]=0;
//  vector_set_used(matches->qbuf_keys_buffer,char_used+key_len+1);
//  return key_used;
//}
//// We _delete_ delta characters in s2, which is longer by delta
//GEM_INLINE void find_min_del(
//    ch_t* s1, ch_t* s2, const idx_t len,
//    const idx_t delta, idx_t* pos, idx_t* misms) {
//  assert(len>=delta);
//  register idx_t mismatches = 0, i;
//  for (i = delta; i < len; ++i) {
//    if (s1[i - delta] != s2[i]) ++mismatches;
//  }
//  register const idx_t top = len - delta;
//  register idx_t min_mismatches = mismatches, min_del_pos = 0;
//  for (i = 0; i < top; ++i) {
//    if (s1[i] != s2[i + delta]) --mismatches;
//
//    if (s1[i] != s2[i]) ++mismatches;
//
//    if (mismatches < min_mismatches) {
//      min_mismatches = mismatches;
//      min_del_pos = i + 1;
//    }
//  }
//  *pos = min_del_pos;
//  *misms = min_mismatches;
//  return;
//}
///* */
//#define MATCH_GET_MISMATCHES_BODY
//  return ((mismatch*)vector_get_mem(matches->rbuf_mismatches))+match->displacement
//#define MATCHES_ITERATOR_NEXT_BODY
//  register pos_match_t* current_match=iterator->current_match;
//  while (iterator->num_misms<iterator->payload->max_decoded_stratum) {
//    while (iterator->position<vector_get_used(iterator->payload->rbuf_pos)) {
//      if (current_match->mismatches==iterator->num_misms) {
//        iterator->current_match=current_match+1;
//        iterator->position++;
//        return current_match;
//      }
//      iterator->position++;
//      current_match++;
//    }
//    iterator->position=0;
//    current_match=vector_get_mem(iterator->payload->rbuf_pos);
//    iterator->num_misms++;
//  }
//  return 0
///* */
//// This one is still private (it provides unsafe non-const access)
//GEM_INLINE pos_match_t* matches_iterator_next_nonconst(matches_iterator* const iterator)
//{ MATCHES_ITERATOR_NEXT_BODY; }
//// This one as well (we do not know about the direction at this level)
//GEM_INLINE direction_t match_get_direction(
//    const pos_match_t* const match,const matches* const mtches)
//{ return ((keys_info*)vector_get_mem(mtches->qbuf_keys_info)+match->key_id)->direction; }
//
////
//// Public functions
////
//GEM_TOPLEVEL_INLINE matches* matches_new() {
//  matches* res=malloc(sizeof(matches));
//  gem_cond_fatal_error(!res,MEM_HANDL);
//  res->rbuf_int=vector_new(FMI_MATCHES_INIT_INTERVALS,sizeof(int_match_t));
//  res->num_int_matches=0;
//  res->rbuf_pos=vector_new(FMI_MATCHES_INIT_POS,sizeof(pos_match_t));
//  res->rbuf_mismatches=vector_new(FMI_MATCHES_INIT_AUG_MISMATCHES*
//                                  FMI_MATCHES_INIT_INTERVALS,sizeof(mismatch));
//  res->rbuf_counts=vector_new(GEM_QUERY_INIT_LEN,sizeof(idx_t));
//  res->max_decoded_stratum=0;
//  res->qbuf_keys_buffer=vector_new(FMI_MATCHES_NUM_QUERIES*GEM_QUERY_INIT_LEN,sizeof(ch_t));
//  res->qbuf_keys_info=vector_new(FMI_MATCHES_NUM_QUERIES,sizeof(keys_info));
//  /* Safe */
//  res->safe_rbuf_counts=vector_new(GEM_QUERY_INIT_LEN,sizeof(idx_t));
//  res->safe_num_int_matches=0;
//  return res;
//}
//GEM_TOPLEVEL_INLINE void matches_delete(matches* const m) {
//  vector_delete(m->rbuf_int);
//  vector_delete(m->rbuf_pos);
//  vector_delete(m->rbuf_mismatches);
//  vector_delete(m->rbuf_counts); //FIXME key buffers
//  free(m);
//  return;
//}
//GEM_INLINE void fm_matches_store_state(matches* const matches) {
//  register const uint64_t num_counts = vector_get_used(matches->rbuf_counts);
//  register const uint64_t* const counts = vector_get_mem(matches->rbuf_counts);
//  vector_reserve(matches->safe_rbuf_counts,num_counts);
//  register uint64_t* const safe_counts = vector_get_mem(matches->safe_rbuf_counts);
//  register uint64_t i;
//  for (i=0; i<num_counts; ++i) {
//    safe_counts[i] = counts[i];
//  }
//  vector_set_used(matches->safe_rbuf_counts,num_counts);
//  matches->safe_rbuf_int_used=vector_get_used(matches->rbuf_int);
//  matches->safe_num_int_matches=matches->num_int_matches;
//  matches->safe_rbuf_pos_used=vector_get_used(matches->rbuf_pos);
//  matches->safe_rbuf_mismatches_used=vector_get_used(matches->rbuf_mismatches);
//  matches->safe_bad_quality_misms=matches->bad_quality_misms;
//  matches->safe_max_complete_stratum=matches->max_complete_stratum;
//}
//GEM_INLINE void fm_matches_recover_state(matches* const matches) {
//  register const uint64_t num_counts = vector_get_used(matches->safe_rbuf_counts);
//  register uint64_t* const counts = vector_get_mem(matches->rbuf_counts);
//  register const uint64_t* const safe_counts = vector_get_mem(matches->safe_rbuf_counts);
//  register uint64_t i;
//  for (i=0; i<num_counts; ++i) {
//    counts[i] = safe_counts[i];
//  }
//  vector_set_used(matches->rbuf_counts,num_counts);
//  vector_set_used(matches->rbuf_int,matches->safe_rbuf_int_used);
//  matches->num_int_matches=matches->safe_num_int_matches;
//  vector_set_used(matches->rbuf_pos,matches->safe_rbuf_pos_used);
//  vector_set_used(matches->rbuf_mismatches,matches->safe_rbuf_mismatches_used);
//  matches->bad_quality_misms=matches->safe_bad_quality_misms;
//  matches->max_complete_stratum=matches->safe_max_complete_stratum;
//}
//GEM_TOPLEVEL_INLINE idx_t matches_get_length(const matches* const mtches) {
//  register const vector* counts=mtches->rbuf_counts;
//  register idx_t* cntrs=vector_get_mem(counts),acc=0,i;
//  for (i=0;i<vector_get_used(counts);++i) acc+=cntrs[i];
//  return acc;
//}
//GEM_TOPLEVEL_INLINE const vector* matches_get_counts(const matches* const mtches)
//{ return mtches->rbuf_counts; }
//GEM_TOPLEVEL_INLINE void matches_clear(matches* const m) {
//  vector_set_used(m->rbuf_int,0);
//  vector_set_used(m->rbuf_pos,0);
//  vector_set_used(m->rbuf_mismatches,0);
//  vector_set_used(m->rbuf_counts,0);
//  vector_set_used(m->qbuf_keys_info,0);
//  vector_set_used(m->qbuf_keys_buffer,0);
//  m->num_int_matches=0;
//  m->bad_quality_misms=0;
//  m->max_decoded_stratum=0;
//  m->max_complete_stratum=UINT64_MAX;
//  /* Safe */
//  vector_set_used(m->safe_rbuf_counts,0);
//  m->safe_rbuf_int_used=0;
//  m->safe_rbuf_pos_used=0;
//  m->safe_rbuf_mismatches_used=0;
//  m->safe_bad_quality_misms=0;
//  m->safe_max_complete_stratum=UINT64_MAX;
//  return;
//}
//GEM_TOPLEVEL_INLINE
//void matches_get_iterator(const matches* const mtches,matches_iterator* iter) {
//  iter->payload=mtches;
//  iter->num_misms=0;
//  iter->position=0;
//  iter->current_match=vector_get_mem(mtches->rbuf_pos);
//  return;
//}
//GEM_TOPLEVEL_INLINE bool matches_iterator_is_empty(const matches_iterator* const iterator)
//{ return iterator->payload->max_decoded_stratum==0; }
//GEM_TOPLEVEL_INLINE const pos_match_t* matches_iterator_next(matches_iterator* const iterator)
//{ MATCHES_ITERATOR_NEXT_BODY; }
//GEM_TOPLEVEL_INLINE uint64_t match_get_position(const pos_match_t* const match)
//{ return match->position; }
//GEM_TOPLEVEL_INLINE uint64_t match_get_query_id(const pos_match_t* const match)
//{ return match->key_id; }
//GEM_TOPLEVEL_INLINE const ch_t* match_get_query(
//    const pos_match_t* const match,const matches* const mtches)
//{ return (ch_t*)vector_get_mem(mtches->qbuf_keys_buffer)+
//         ((keys_info*)vector_get_mem(mtches->qbuf_keys_info)+match->key_id)->displacement; }
//GEM_TOPLEVEL_INLINE uint64_t match_get_query_length(
//    const pos_match_t* const match,const matches* const mtches)
//{ return ((keys_info*)vector_get_mem(mtches->qbuf_keys_info)+match->key_id)->key_len; }
//GEM_TOPLEVEL_INLINE uint64_t match_get_mismatches_length(const pos_match_t* const match) {
//  return match->mismatches;
//}
//GEM_TOPLEVEL_INLINE const mismatch* match_get_mismatches( // Public const version
//    const pos_match_t* const match,const matches* const matches) {
//  MATCH_GET_MISMATCHES_BODY;
//}
//GEM_INLINE mismatch* match_get_mismatches_nonconst( // Private non-const version
//    const pos_match_t* const match,const matches* const matches) {
//  MATCH_GET_MISMATCHES_BODY;
//}
///* */
//GEM_TOPLEVEL_INLINE multimatches* multimatches_new() {
//  multimatches* mmatches = malloc(sizeof(multimatches));
//  gem_cond_fatal_error(!mmatches,MEM_HANDL);
//  mmatches->matches_groups=vector_new(FMI_MULTIMATCHES_INIT_GROUPS,sizeof(matches*));
//  mmatches->multimatches=vector_new(FMI_MULTIMATCHES_INIT_MATCHES,sizeof(multimatch));
//  mmatches->mm_counts=vector_new(GEM_QUERY_INIT_LEN,sizeof(idx_t));
//  return mmatches;
//}
//GEM_TOPLEVEL_INLINE void multimatches_delete(multimatches* const mmatches) {
//  vector_delete(mmatches->matches_groups);
//  vector_delete(mmatches->multimatches);
//  vector_delete(mmatches->mm_counts);
//  free(mmatches);
//}
//GEM_TOPLEVEL_INLINE void multimatches_clear(multimatches* const mmatches) {
//  vector_set_used(mmatches->multimatches,0);
//  vector_set_used(mmatches->mm_counts,0);
//}
//GEM_TOPLEVEL_INLINE void multimatches_add_group(
//    multimatches* const multimatches,const matches* const mtchs) {
//  vector_reserve(multimatches->matches_groups,vector_get_used(multimatches->matches_groups));
//  register const matches** ptr_matches = vector_get_mem(multimatches->matches_groups);
//  ptr_matches[vector_get_used(multimatches->matches_groups)] = mtchs;
//  vector_inc_used(multimatches->matches_groups);
//}
//
//
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//#define FMI_MATCHES_ADJUST_MAX_MISMS_USING_DELTA(distance)
//  if (search_params->max_mismatches > search_params->delta+distance) {
//    search_params->max_mismatches = search_params->delta+distance;
//  }
//
//#define FMI_MATCHES_STORE_MATCH(pos_match,match_position,match_length,match_num_misms)
//  pos_match->position = match_position;
//  pos_match->length = match_length;
//  pos_match->mismatches = match_num_misms;
//  pos_match->displacement = misms - ((mismatch*) vector_get_mem(matches->rbuf_mismatches));
//  pos_match->key_id = search_params->key_id;
//  FMI_MATCHES_ADJUST_MAX_MISMS_USING_DELTA(match_num_misms)
//
//
//#define FMI_MATCHES_STORE_NEW_MATCH(pos_match,position,length,num_misms)
//  FMI_MATCHES_STORE_MATCH(pos_match,position,length,num_misms);
//  misms+=num_misms;
//  ++cntrs[num_misms];
//  ++pos_match
//
//#define FMI_MATCHES_REPLACE_MATCH(pos_match,position,length,old_num_misms,new_num_misms)
//  --cntrs[old_num_misms];
//  ++cntrs[new_num_misms];
//  FMI_MATCHES_STORE_MATCH(pos_match,position,length,new_num_misms);
//  misms+=new_num_misms
//
//#define FMI_MATCHES_CHECK_DUPLICATES__STORE_MATCH(match_position,match_length,match_num_misms) {
//  register pos_match_t* pos_match_it = (pos_match_t*)vector_get_mem(matches->rbuf_pos);
//  register const uint64_t num_pos_matches = out-pos_match_it;
//  register bool is_duplicated;
//  register uint64_t it;
//  for (it=0, is_duplicated=false; !is_duplicated && it<num_pos_matches; ++it, ++pos_match_it) {
//    if (MATCHES_IS_DUPLICATED(pos_match_it->position,match_position,
//        pos_match_it->length,match_length,pos_match_it->key_id,search_params->key_id)) {
//      is_duplicated = true;
//      if (pos_match_it->mismatches < match_num_misms) break;
//      if (pos_match_it->mismatches == match_num_misms) {
//        register mismatch* misms_current_match = (mismatch*)vector_get_mem(matches->rbuf_mismatches)+pos_match_it->displacement;
//        register uint64_t lev_current_match=0, misms_it;
//        register uint64_t lev_new_match=0;
//        for (misms_it=0;misms_it<match_num_misms;++misms_it) {
//          lev_current_match += FMI_MATCHES_IS_INDEL(misms_current_match[misms_it].mismatch) ?
//              FMI_MATCHES_GET_SIZE_INDEL(misms_current_match[misms_it].mismatch) : 1;
//          lev_new_match += FMI_MATCHES_IS_INDEL(misms[misms_it].mismatch) ?
//              FMI_MATCHES_GET_SIZE_INDEL(misms[misms_it].mismatch) : 1;
//        }
//        if (lev_current_match <= lev_new_match) break;
//      }
//      FMI_MATCHES_REPLACE_MATCH(pos_match_it,match_position,match_length,
//          pos_match_it->mismatches,match_num_misms);
//    }
//  }
//  if (!is_duplicated) {
//    FMI_MATCHES_STORE_NEW_MATCH(out,match_position,match_length,match_num_misms);
//  }
//}
//
//// Myers DEBUG
//bool fmi_matches_check_alignment(
//    const ch_t* const key,const uint64_t key_length,
//    const ch_t* const text,uint64_t alg_length,
//    mismatch* misms,const int64_t num_misms) {
//  register bool correct = true;
//  register int64_t p_ref = 0, p_read = 0, misms_pos=0;
//  while (correct && p_read<key_length && p_ref<alg_length) {
//    if (misms_pos<num_misms && misms[misms_pos].position == p_read) {
//      if (misms[misms_pos].mismatch>=256) {
//        register const uint64_t tmp=misms[misms_pos].mismatch/256;
//        register const bool is_del = tmp%2;
//        register const uint64_t size = tmp/2;
//        if (!is_del) p_ref+=size;
//        else p_read+=size;
//      } else {
//        correct = (misms[misms_pos].mismatch == text[p_ref] && key[p_read] != text[p_ref]);
//        p_read++; p_ref++;
//      }
//      ++misms_pos;
//    } else {
//      correct = key[p_read] == text[p_ref];
//      p_read++; p_ref++;
//    }
//  }
//  while (correct && p_read<key_length) {
//    if (misms_pos>=0 && misms[misms_pos].position == p_read) {
//      if (misms[misms_pos].mismatch>=256) {
//        register const uint64_t tmp=misms[misms_pos].mismatch/256;
//        register const bool is_del = tmp%2;
//        register const uint64_t size = tmp/2;
//        if (!is_del) correct = false;
//        else p_read+=size;
//      } else {
//        correct = false;
//      }
//      ++misms_pos;
//    } else {
//      correct = key[p_read] == text[p_ref];
//      p_read++; p_ref++;
//    }
//  }
//  return correct && p_read==key_length && misms_pos==num_misms && p_ref==alg_length;
//}
//
//GEM_TOPLEVEL_INLINE bool fmi_matches_check_levenshtein_match(
//    const ch_t* const key,const uint64_t* const peq,const uint64_t key_length,
//    const ch_t* const decoded_text,const uint64_t text_length,
//    const uint64_t distance,const bool allow_indels,
//    mismatch* const misms,idx_t* match_position,idx_t* match_end_position,
//    uint64_t* match_distance,vector_pool* const mpool) {
//  register const uint64_t key_mod = key_length%WORD_SIZE_64;
//  register const uint64_t num_words = (key_length+(WORD_SIZE_64-1))/WORD_SIZE_64;
//  register const uint64_t scope = FMI_MATCHES_COMMON_SCOPE;
//  int64_t* score;
//  uint64_t* level_mask;
//  int64_t* init_score;
//  register uint64_t i, j;
//  uint64_t* vP;
//  uint64_t* vM;
//
//  // Allocate memory
//  FMI_MATCHES_DP_PREPARE_MEM_M();
//
//  // Init DP structures
//  register uint8_t y;
//  FMI_MATCHES_DP_INIT_M(0);
//
//  register uint64_t hit_pos, score_match, opt_positions;
//  register bool match_found=false;
//  register int8_t carry;
//  for (j=0; j<text_length; ++j) {
//    register const uint8_t enc_char = bwt_dna_p_encode[decoded_text[j]];
//    register const uint64_t current_pos = j;
//    register const uint64_t next_pos = j+1;
//    // Advance blocks and check cut off strategy
//    FMI_MATCHES_DP_ADVANCE();
//    FMI_MATCHES_DP_CUT_OFF();
//    // Check match and its optimization
//    FMI_MATCHES_CHECK__OPT_MATCH(text_length-1);
//  }
//
//  // Retrieve the alignment
//  if (match_found) {
//    register bool has_indel = false;
//    // Store the match (and backtrace the mismatches)
//    FMI_MATCHES_BACKTRACE__STORE_MATCH(0,h+1,h);
//    // if (!allow_indels && !has_indel) return false; // Avoid hamming duplicates
//    *match_position = h+1;
//    *match_end_position = hit_pos+1;
//    *match_distance = score_match-misms_num;
//    ADD_COUNTER(GSC_ONLINE_ASM_MISMS, *match_distance);
//    ADD_COUNTER(GSC_ONLINE_ASM_ERRORS, score_match);
//#ifdef GEM_MAPPER_DEBUG // Check alignment // FIXME: DON'T
//    if (!fmi_matches_check_alignment(key,key_length,decoded_text+(*match_position),hit_pos-(*match_position)+1,misms,score_match-misms_num)) {
//      //fprintf(stderr,"Wrong alignment\n"); exit(0);
//      return false;
//    }
//#endif
//    return true;
//  }
//  return false;
//}
//
//#define FMI_MATCHES_DP_INIT_V()
//  y = (distance>0) ? (distance+(WORD_SIZE_64-1))/WORD_SIZE_64 : 1;
//  score[0] = init_score[0];
//  for (i=1; i<y; ++i) score[i] = score[i-1] + init_score[i];
//  for (i=0; i<y; ++i) INIT_BLOCK_V_(i)
//#define FMI_MATCHES_SEARCH_NEXT_CHAR(enc_char,forward_search)
//  if (forward_search) {
//    enc_char = be_stream_char_encoded(be,&be_stream);
//    decoded_text[j%scope] = bwt_dna_p_decode[enc_char];
//  } else {
//    enc_char = be_inverse_stream_char_encoded(be,&be_stream);
//    decoded_text[(sum_limits-j)%scope] = bwt_dna_p_decode[enc_char];
//  }
//GEM_TOPLEVEL_INLINE uint64_t fmi_matches_search(
//    const _FMI_* const index,fmi_extend_parameters* const extend_parameters,
//    const uint64_t init_pos, const uint64_t final_pos,
//    const int64_t distance, const uint64_t max_num_matches,
//    const bool* const allowed_repl,matches* const matches,vector_pool* const mpool) {
//  // Pattern vars
//  register const uint64_t key_length = extend_parameters->key_len;
//  register const bool forward_search = (extend_parameters->search_direction==Forward);
//  register const uint64_t* peq = forward_search ?
//      extend_parameters->forward_peq : extend_parameters->reverse_peq;
//
//  // Search vars
//  register const uint64_t key_mod = key_length%WORD_SIZE_64;
//  register const uint64_t num_words = (key_length+(WORD_SIZE_64-1))/WORD_SIZE_64;
//  register const uint64_t scope = FMI_MATCHES_COMMON_SCOPE;
//  register uint64_t i, j, num_matches=0;
//
//  // Matches vars
//  register const uint64_t possible_matches = UNSAFE_MIN((final_pos-init_pos+1+key_length)/key_length,max_num_matches);
//  matches_prepare(matches,distance);
//  register uint64_t* cntrs = vector_get_mem(matches->rbuf_counts);
//  vector_reserve(matches->rbuf_pos,vector_get_used(matches->rbuf_pos)+possible_matches);
//  register pos_match_t* pos_match = vector_get_mem_next_elm(matches->rbuf_pos,pos_match_t);
//  vector_reserve(matches->rbuf_mismatches,vector_get_used(matches->rbuf_mismatches)+distance*possible_matches);
//  register mismatch* misms = vector_get_mem_next_elm(matches->rbuf_mismatches,mismatch);
//  idx_t match_position, match_end_position;
//  uint64_t match_distance;
//
//  // Allocate memory
//  vector_prepare(mpool->buf1,uint64_t,5*num_words);
//  uint64_t* P = vector_get_mem(mpool->buf1);
//  uint64_t* M = P + num_words;
//  uint64_t* level_mask = M + num_words;
//  int64_t* score = (int64_t*)level_mask + num_words;
//  int64_t* init_score = score + num_words;
//  vector_prepare(mpool->fbuf2,ch_t,2*scope);
//  ch_t* decoded_text = vector_get_mem(mpool->fbuf2);
//  ch_t* match_text = decoded_text+scope;
//
//  // Init
//  register const uint8_t top = num_words-1;
//  for (i=0; i<top; ++i) {
//    level_mask[i] = WORD_MASK_64;
//    init_score[i] = WORD_SIZE_64;
//  }
//  level_mask[top] = (key_mod>0) ? 1L<<(key_mod-1) : WORD_MASK_64;
//  init_score[top] = (key_mod>0) ? key_mod : WORD_SIZE_64;
//  register uint8_t y;
//  FMI_MATCHES_DP_INIT_V();
//
//  // Advance in DP-bit_encoded matrix
//  register const BE_DNA* const be = index->bed;
//  register int8_t carry;
//  BE_stream_DNA be_stream;
//  if (forward_search) {
//    be_stream_new(be,init_pos,&be_stream);
//  } else {
//    be_inverse_stream_new(be,final_pos,&be_stream);
//  }
//  register const uint64_t sum_limits = init_pos+final_pos;
//  for (j=init_pos; j<=final_pos; ++j) {
//    // Fetch next character
//    register slch_t enc_char;
//    FMI_MATCHES_SEARCH_NEXT_CHAR(enc_char,forward_search);
//    // Reset if not allowed character
//    if (__builtin_expect(!allowed_repl[enc_char],false)) {
//      if (__builtin_expect(enc_char==CHAR_ENC_SEP||enc_char==CHAR_ENC_EOT,false)) return num_matches;
//      FMI_MATCHES_DP_INIT_V(); continue;
//    }
//    // Advance all blocks
//    for (i=0,carry=0; i<y; ++i) {
//      register uint64_t* const Py = P+i;
//      register uint64_t* const My = M+i;
//      carry = advance_block(peq[VDP_IDX(enc_char,i,num_words)],level_mask[i],*Py,*My,carry+1,Py,My);
//      score[i] += carry;
//    }
//
//    // Cut-off
//    register const uint8_t last = y-1;
//    if ((score[last]-carry) <= distance && last<top &&
//        ( (peq[PEQ_IDX(enc_char,y,num_words)] & 1) || (carry<0) )  ) {
//      INIT_BLOCK_V_(y);
//      register uint64_t* const Py = P+y;
//      register uint64_t* const My = M+y;
//      score[y] = score[y-1] + init_score[y] - carry +
//          advance_block(peq[VDP_IDX(enc_char,y,num_words)],level_mask[y],*Py,*My,carry+1,Py,My);
//      ++y;
//    } else {
//      while (score[y-1] > distance+init_score[y-1]) {
//        --y;
//      }
//    }
//
//    // Check match
//    if (y==num_words && score[y-1]<=distance) {
//      INC_COUNTER(GSC_ONLINE_ASM_HIT);
//      ADD_COUNTER(GSC_ONLINE_ASM_DISTANCE,(j-init_pos+key_length));
//      register const uint64_t offset = scope-1;
//      register uint64_t start_pos, end_pos;
//      // Extend @distance characters to optimize the alignment
//      register const uint64_t extended_limit = j+distance<=final_pos ? j+distance : final_pos;
//      for (++j;j<=extended_limit;++j) {
//        register slch_t enc_char;
//        FMI_MATCHES_SEARCH_NEXT_CHAR(enc_char,forward_search);
//        if (__builtin_expect(!allowed_repl[enc_char],false)) break;
//      }
//      --j;
//      if (forward_search) { // Interval alignment
//        start_pos = (j-init_pos)>=offset?j-offset:init_pos;
//        end_pos = j;
//      } else {
//        start_pos = (sum_limits-j);
//        end_pos = (start_pos+offset)<=final_pos?start_pos+offset:final_pos;
//      }
//      // Copy the decoded text and extract the alignment
//      for (i=start_pos; i<=end_pos; ++i) {
//        match_text[i-start_pos]=decoded_text[i%scope];
//      }
////      if(!fmi_matches_retrieve_dp_match(index,extend_parameters->key_id,
////          extend_parameters->forward_key,extend_parameters->forward_peq,
////          key_length,start_pos,end_pos,decoded_text,distance,matches,mpool)) { // FIXME: DEBUG
////        printf("Bad STD align\n");
////      }
//      gem_cond_fatal_error(!fmi_matches_check_levenshtein_match(
//          extend_parameters->forward_key,extend_parameters->forward_peq,key_length,
//          match_text,end_pos-start_pos+1,distance,true,
//          misms,&match_position,&match_end_position,&match_distance,mpool),
//          ALGO_REPORT,"Bad paired alignment");
//      pos_match->position = match_position+start_pos;
//      pos_match->length = match_end_position-match_position;
//      pos_match->mismatches = match_distance;
//      pos_match->displacement = misms - ((mismatch*) ((matches->rbuf_mismatches)->memory));
//      pos_match->key_id = extend_parameters->key_id;
//      misms+=match_distance;
//      ++cntrs[match_distance]; ++pos_match; ++num_matches;
//      vector_update_used(matches->rbuf_mismatches,misms);
//      vector_update_used(matches->rbuf_pos,pos_match);
//      if (num_matches==max_num_matches) return num_matches;
//      FMI_MATCHES_DP_INIT_V();
//    }
//  }
//  return num_matches;
//}
//
//inline bool fmi_check_single_indel_candidate(
//    const _FMI_* const a, ch_t* key,const idx_t pos,
//    ch_t* decoded_key,const uint64_t key_len,const uint64_t delta,
//    const bool* const allowed_chars,const ch_t* mismatch_mask,
//    const uint64_t max_mismatches,idx_t *effective_position,
//    idx_t *effective_end_position,uint64_t *total_mismatches,mismatch* misms_info) {
//  register idx_t j;
//  register bool is_del = true;
//  idx_t mism, new_mism, gap_pos, new_gap_pos;
//  register idx_t eff_pos = pos - delta;
//  register mismatch* misms_info_start = misms_info;
//
//  register ch_t* indeled_decoded = decoded_key - delta;
//
//  find_min_del(indeled_decoded, key, key_len, delta, &gap_pos, &mism);
//  find_min_del(key, indeled_decoded, key_len + delta, delta, &new_gap_pos, &new_mism);
//  if (new_mism < mism) {
//    is_del = false;
//    mism = new_mism;
//    gap_pos = new_gap_pos;
//  }
//  find_min_del(decoded_key, key, key_len, delta, &new_gap_pos, &new_mism);
//  if (new_mism < mism) {
//    is_del = true;
//    mism = new_mism;
//    gap_pos = new_gap_pos;
//    eff_pos = pos;
//    indeled_decoded = decoded_key;
//  }
//  find_min_del(key, decoded_key, key_len + delta, delta, &new_gap_pos,&new_mism);
//  if (new_mism < mism) {
//    is_del = false;
//    gap_pos = new_gap_pos;
//    eff_pos = pos;
//    indeled_decoded = decoded_key;
//  }
//
//  // Check if the indel is equivalent to a mismatch or is it a real indel
//  if (delta >= 1 && (gap_pos == 0 || gap_pos == key_len - 1)) return false;
//
//  // We do not know the number of mismatches yet
//  register idx_t num_mismatches = 0;
//  for (j = 0; j < gap_pos; ++j) {
//    if (indeled_decoded[j] != key[j]) {
//      if (!allowed_chars[indeled_decoded[j]]) return false;
//      if (allowed_chars[key[j]] && (mismatch_mask==0 || mismatch_mask[j]==Real)) {
//        if (++num_mismatches > max_mismatches) return false;
//      }
//
//      misms_info->position = j;
//      misms_info->mismatch = indeled_decoded[j];
//      ++misms_info;
//    }
//  }
//
//  /* The indel itself.
//     We do not increment mismatches here, since we consider the indel as
//     one of the additional mismatches allowed for bad-quality bases, and
//     we do not want it to contribute to the selection; however, it will
//     count when it is inserted later on */
//  misms_info->position = gap_pos;
//  misms_info->mismatch = FMI_MATCHES_INDEL_VALUE(is_del,delta);
//  ++misms_info;
//
//  // Check the rest of the decoded read (after the indel)
//  key += gap_pos;
//  if (mismatch_mask) mismatch_mask += gap_pos;
//  indeled_decoded += gap_pos;
//  register idx_t top = key_len - gap_pos;
//
//  if (is_del) {
//    key += delta;
//    if (mismatch_mask) mismatch_mask += delta;
//    top -= delta;
//  } else {
//    indeled_decoded += delta;
//  }
//
//  for (j = 0; j < top; ++j) {
//    if (indeled_decoded[j] != key[j]) {
//      if (!allowed_chars[indeled_decoded[j]]) return false;
//      if (allowed_chars[key[j]] && (mismatch_mask==0 || mismatch_mask[j]==Real)) {
//        if (++num_mismatches > max_mismatches) return false;
//      }
//
//      misms_info->position = j + gap_pos + (is_del ? delta : 0);
//      misms_info->mismatch = indeled_decoded[j];
//      ++misms_info;
//    }
//  }
//
//  // All the mismatches
//  *total_mismatches = misms_info - misms_info_start;
//  *effective_position = eff_pos;
//  *effective_end_position = eff_pos + key_len + (is_del ? -delta : +delta);
//  return true;
//}
//
///*
// * Batch decode of all positions of the candidates => end_key_pos = decoded(end_reg_idx);
// * ** CHECKED[12/7/2011]
// */
//GEM_INLINE
//void fmi_matches_decode_positions_candidates(const _FMI_* const a,const uint64_t key_len,
//                                             matches* const matches,
//                                             const vector* const positions) {
//  register uint64_t i;
//  register idx_t pos;
//  register filter_query* queries = vector_get_mem(positions);
//  for (i=0; i<vector_get_used(positions); ++i) {
//    pos = fmi_lookup(a, queries[i].end_reg_idx);
//    if (pos >= queries[i].end) {
//      pos -= queries[i].end;
//      if ((pos + key_len) < a->bwt->n) { // FIXME: CheckMe a->bwt->n? or text_length
//        queries[i].end_key_pos = pos;
//      } else {
//        queries[i].end_key_pos = IDX_T_MAX;
//      }
//    } else {
//      queries[i].end_key_pos = IDX_T_MAX;
//    }
//  }
//}
//
//#define GET_NEXT_LOOKUP_QUERY(placeholder)
//  next_pos[placeholder]=queries[count].end_reg_idx;
//  id[placeholder]=count;
//  ++count;
//  dist[placeholder]=0;
//  is_sampled[placeholder]=false
//
///*
// * Bach decode of all filter queries and aligned to the beginning of the read.
// * All the steps (CSA-lookup, rankQueries) are performed with prefetch-loops
// * ** CHECKED[12/7/2011]
// */
//GEM_INLINE
//void fmi_matches_prefetched_decode_positions_candidates(
//    const _FMI_* const a,const uint64_t key_len,
//    const matches* const matches,vector* const positions) {
//  register filter_query* queries = vector_get_mem(positions);
//  register uint64_t count, num_left, i;
//
//  idx_t id[NUMBER_DECODE_QUERIES_PREFETCHED];
//  idx_t dist[NUMBER_DECODE_QUERIES_PREFETCHED];
//  idx_t pos[NUMBER_DECODE_QUERIES_PREFETCHED];
//  idx_t next_pos[NUMBER_DECODE_QUERIES_PREFETCHED];
//  bool is_sampled[NUMBER_DECODE_QUERIES_PREFETCHED];
//  bwt_prefetch_info call_info[NUMBER_DECODE_QUERIES_PREFETCHED];
//
//  // Sampled-LF loop. Fill the buffer with queries
//  num_left = UNSAFE_MIN(vector_get_used(positions), NUMBER_DECODE_QUERIES_PREFETCHED);
//  for (i=0, count=0; i<num_left; ++i) {
//    GET_NEXT_LOOKUP_QUERY(i);
//  }
//  for (; i<NUMBER_DECODE_QUERIES_PREFETCHED; ++i) {
//    is_sampled[i] = true;
//  }
//  // Full-prefetch loop for sampled-LF
//  if (num_left == NUMBER_DECODE_QUERIES_PREFETCHED) {
//    while (count < vector_get_used(positions)) {
//      for (i=0; i<NUMBER_DECODE_QUERIES_PREFETCHED; ++i) {
//        pos[i] = next_pos[i];
//        bwt_sampled_LF_prefetch(a->bwt,pos[i],call_info+i);
//      }
//      for (i=0; i<NUMBER_DECODE_QUERIES_PREFETCHED; ++i) {
//        next_pos[i] = bwt_sampled_LF_bit__erank_retrieve(a->bwt,pos[i],is_sampled+i,call_info+i);
//        ++dist[i];
//        if (is_sampled[i]) {
//          queries[id[i]].f_sampled_pos = pos[i];
//          queries[id[i]].f_distance = dist[i];
//          if (count < vector_get_used(positions)) {
//            GET_NEXT_LOOKUP_QUERY(i);
//          }
//        }
//      }
//    }
//  }
//  // Solve remaining queries
//  for (i=0; i<NUMBER_DECODE_QUERIES_PREFETCHED; ++i) {
//    if (is_sampled[i]) continue;
//    do {
//      pos[i] = next_pos[i];
//      next_pos[i] = bwt_sampled_LF_bit__erank(a->bwt,pos[i],is_sampled+i);
//      ++dist[i];
//    } while (!is_sampled[i]);
//    queries[id[i]].f_sampled_pos = pos[i];
//    queries[id[i]].f_distance = dist[i];
//  }
//
//  // Prefetch LF-erank
//  count = vector_get_used(positions); // Indeed this is already the current value
//  num_left = NUMBER_DECODE_QUERIES_PREFETCHED;
//  while (num_left == NUMBER_DECODE_QUERIES_PREFETCHED) {
//    num_left = UNSAFE_MIN(count, NUMBER_DECODE_QUERIES_PREFETCHED);
//    for (i=0; i<num_left; ++i) {
//      bwt_sampled_LF_prefetch(a->bwt,queries[i].f_sampled_pos,call_info+i);
//    }
//    for (i=0; i<num_left; ++i, ++queries) {
//      queries->f_csa_pos = bwt_sampled_LF_erank_retrieve(
//          a->bwt,queries->f_sampled_pos,call_info+i);
//    }
//    count -= NUMBER_DECODE_QUERIES_PREFETCHED;
//  }
//
//  // Prefetch CSA query
//  uint64_t* addr[NUMBER_DECODE_QUERIES_PREFETCHED];
//  count = vector_get_used(positions);
//  num_left = NUMBER_DECODE_QUERIES_PREFETCHED;
//  queries = vector_get_mem(positions);
//  while (num_left == NUMBER_DECODE_QUERIES_PREFETCHED) {
//    num_left = UNSAFE_MIN(count, NUMBER_DECODE_QUERIES_PREFETCHED);
//    for (i=0; i<num_left; ++i) {
//      addr[i] = GEM_SA_ADDRESS(a->D,a->cntr_bytes,queries[i].f_csa_pos);
//      GEM_PREFETCH(addr[i]);
//    }
//    for (i=0; i<num_left; ++i, ++queries) {
//      queries->end_key_pos = (GEM_UNALIGNED_ACCESS_GIVEN_ADDR(addr[i],a->cntr_bytes)
//          + queries->f_distance - 1) % (a->text_length+1);
//    }
//    count -= NUMBER_DECODE_QUERIES_PREFETCHED;
//  }
//
//  // Adjust decoded position to the beginning of the read
//  queries = vector_get_mem(positions);
//  for (i=0; i<vector_get_used(positions); ++i, ++queries) {
//    if (queries->end_key_pos >= queries->end) {
//      queries->end_key_pos -= queries->end;
//      if ((queries->end_key_pos + key_len) >= a->bwt->n) {
//        queries->end_key_pos = IDX_T_MAX;
//      }
//    } else {
//      queries->end_key_pos = IDX_T_MAX;
//    }
//  }
//}
//#define CANDIDATE_CHECK_CHARACTER(CHAR,IDX,MISMS,MISMS_NUM)
//  if (CHAR != key[IDX]) {
//    if (!allowed_chars[CHAR]) return false;
//    if (allowed_chars[key[IDX]] && (mismatch_mask==0 || mismatch_mask[IDX]==Real)) {
//      if (++MISMS_NUM > max_mismatches) return false;
//      misms_ratio+=2;
//      if (use_density && misms_ratio>max_mismatches) return false;
//    } else if (misms_ratio > 0) {
//      --misms_ratio;
//    }
//    MISMS->position = IDX;
//    MISMS->mismatch = CHAR;
//    ++MISMS;
//  }
//GEM_TOPLEVEL_INLINE bool fmi_matches_check_hamming_match(
//    const ch_t* const key,const uint64_t key_length,
//    const ch_t* const mismatch_mask,const bool* const allowed_chars,
//    const ch_t* const decoded_text,const uint64_t max_mismatches,
//    mismatch* misms,uint64_t* const match_distance,vector_pool* const mpool) {
//  register uint64_t i, mismatches, misms_ratio;
//  register mismatch* misms_init = misms;
//  register const bool use_density = max_mismatches>3;
//  for (i=0,mismatches=0,misms_ratio=0; i<key_length; ++i) {
//    CANDIDATE_CHECK_CHARACTER(decoded_text[i],i,misms,mismatches);
//  }
//  *match_distance = misms - misms_init;
//  return true;
//}
//
//#define FMI_MATCHES_ADD_INITIAL_TRIM(misms,length_trim,match_position,match_mismatches) {
//  register int64_t j;
//  if (match_mismatches>0 && misms[0].position == 0 &&
//      FMI_MATCHES_IS_INDEL(misms[0].mismatch) &&
//      FMI_MATCHES_IS_DEL(misms[0].mismatch)) {
//    misms[0].mismatch = FMI_MATCHES_INDEL_VALUE(true,
//        length_trim+FMI_MATCHES_GET_SIZE_INDEL(misms[0].mismatch));
//    for (j=1; j<match_mismatches; ++j) {
//      misms[j].position += length_trim;
//    }
//  } else {
//    for (j=match_mismatches-1; j>=0; --j) {
//      misms[j+1].mismatch = misms[j].mismatch;
//      misms[j+1].position = misms[j].position+length_trim;
//    }
//    misms[0].position = 0;
//    misms[0].mismatch = FMI_MATCHES_INDEL_VALUE(true,length_trim);
//    ++match_mismatches;
//  }
//}
//#define FMI_MATCHES_ADD_FINAL_TRIM(misms,length_trim,key_len,match_end_position,match_mismatches)
//  register const uint64_t last_misms = match_mismatches-1;
//  match_end_position+=length_trim;
//  if (match_mismatches>0 && FMI_MATCHES_IS_INDEL(misms[last_misms].mismatch) &&
//      FMI_MATCHES_IS_DEL(misms[last_misms].mismatch) &&
//      (misms[last_misms].position+FMI_MATCHES_GET_SIZE_INDEL(misms[last_misms].mismatch)) == key_len) {
//    misms[last_misms].mismatch = FMI_MATCHES_INDEL_VALUE(true,
//        length_trim+FMI_MATCHES_GET_SIZE_INDEL(misms[last_misms].mismatch));
//  } else {
//    misms[match_mismatches].position = key_len;
//    misms[match_mismatches].mismatch = FMI_MATCHES_INDEL_VALUE(true,length_trim);
//    ++match_mismatches;
//  }
//
//
//
////
//// Public functions
////
//#define FMI_MATCHES_STORE_DECODED_INT_POSITION(pos_match,dec_position,dec_length,dec_displacement,dec_num_misms,dec_id_key)
//  pos_match->length=dec_length;
//  pos_match->position=dec_position;
//  pos_match->displacement=dec_displacement;
//  pos_match->mismatches=dec_num_misms; /* All the mismatches */
//  pos_match->key_id=dec_id_key;
//  ++pos_match
//GEM_TOPLEVEL_INLINE void fmi_matches_decode(
//    const _FMI_* const index,const uint64_t num_decoded_strata,
//    const uint64_t first_stratum_threshold,const uint64_t num_decoded_matches,
//    matches* const matches,vector_pool* const mpool) {
//  // First, a trivial check
//  if (num_decoded_strata==0 && num_decoded_matches==0) return;
//
//  // Shrink the counter until the last non-zero stratum
//  register idx_t* cntrs=vector_get_mem(matches->rbuf_counts);
//  register int64_t i = vector_get_used(matches->rbuf_counts)-1;
//  while (i>=0 && cntrs[i]==0) --i;
//  vector_get_used(matches->rbuf_counts) = ++i;
//  if (i==0) return;
//
//  // Calculate the maximum reachable stratum w.r.t counters
//  register uint64_t max_stratum = i;
//
//  // Calculate the first non-zero stratum
//  i = 0;
//  while (i<max_stratum && cntrs[i]==0) ++i;
//  register uint64_t min_nz_stratum = i+1;
//
//  // Count matches and strata and calculate the maximum stratum to decode
//  register int64_t acc;
//  for (i=0,acc=0;i<max_stratum;++i) {
//    acc+=cntrs[i];
//    if (acc>num_decoded_matches) break;
//  }
//  register const uint64_t max_stratum_by_matches = i;
//
//  // Check that we need to decode sth
//  if (num_decoded_strata==0 && max_stratum_by_matches<min_nz_stratum) return;
//
//  // Calculate the maximum stratum to decode by stratum // FIXME: Two mandatory strata up to the limit
//  register uint64_t max_stratum_by_stratum;
//  if (num_decoded_strata==0) {
//    max_stratum_by_stratum = 0;
//  } else if (num_decoded_strata==UINT64_MAX) {
//    max_stratum_by_stratum = max_stratum;
//  } else {
//    max_stratum_by_stratum = min_nz_stratum+num_decoded_strata-1;
//    if (cntrs[min_nz_stratum-1]>first_stratum_threshold) {
//      max_stratum_by_stratum = 0;
//    }
//  }
//
//  // Calculate the maximum stratum to decode
//  register const uint64_t max_decoded_stratum = UNSAFE_MAX(max_stratum_by_matches,max_stratum_by_stratum);
//  matches->max_decoded_stratum = max_decoded_stratum;
//  if (max_decoded_stratum==0) return;
//
//  // Decode the intervals having the correct number of mismatches
//  // Append their content to the list of matches in position space.
//
//  // Count number of matches to decode
//  register int_match_t* interval = vector_get_mem(matches->rbuf_int);
//  for (i=0,acc=0;i<vector_get_used(matches->rbuf_int);++i) {
//    if (interval->mismatches < max_decoded_stratum) {
//      acc+=(interval->hi-interval->lo);
//    }
//    ++interval;
//  }
//  if (acc==0) return;
//
//  // Allocate memory for decoded matches
//  vector_reserve(matches->rbuf_pos,vector_get_used(matches->rbuf_pos)+acc);
//  register pos_match_t* out = (pos_match_t*)vector_get_mem(matches->rbuf_pos)+vector_get_used(matches->rbuf_pos);
//  vector_reserve(matches->rbuf_mismatches,vector_get_used(matches->rbuf_mismatches)+(max_decoded_stratum-1)*acc);
//  register mismatch* misms = (mismatch*)vector_get_mem(matches->rbuf_mismatches)+vector_get_used(matches->rbuf_mismatches);
//  // Decode each occurrence and store the mismatches
//  vector_prepare(mpool->fbuf2, ch_t, GEM_QUERY_INIT_LEN+1);
//  interval = vector_get_mem(matches->rbuf_int);
//  for (i=0;i<vector_get_used(matches->rbuf_int);++i,++interval) {
//    if (interval->mismatches<max_decoded_stratum && interval->lo<interval->hi) {
//      // Gather the key information and prepare buffers
//      register const keys_info* key_info =
//          (keys_info*)vector_get_mem(matches->qbuf_keys_info)+interval->key_id;
//      register const uint64_t key_len = key_info->key_len;
//      register const ch_t* key =
//          (ch_t*)vector_get_mem(matches->qbuf_keys_buffer) + key_info->displacement;
//      vector_reserve(mpool->fbuf2, key_len+1);
//      register ch_t* decoded = vector_get_mem(mpool->fbuf2);
//      // Calculate the misms using the first (the same for all the interval)
//      register mismatch* old_misms = misms;
//      register const uint64_t displacement = misms-((mismatch*)vector_get_mem(matches->rbuf_mismatches));
//      register const idx_t position = fmi_lookup(index,interval->lo);
//      fmi_decode(index,position,key_len,decoded);
//      register int64_t j;
//      for (j=0; j<key_len; ++j) {
//        if (decoded[j] != key[j]) {
//          misms->position=j;
//          misms->mismatch=decoded[j];
//          ++misms;
//        }
//      }
//      FMI_MATCHES_STORE_DECODED_INT_POSITION(out,position,key_len,
//          displacement,interval->mismatches,interval->key_id);
//      assert((misms-old_misms)==interval->mismatches);
//      // Add the rest decoded interval positions
//      for (j=interval->lo+1;j<interval->hi;++j) {
//        FMI_MATCHES_STORE_DECODED_INT_POSITION(out,fmi_lookup(index,j),key_len,
//            displacement,interval->mismatches,interval->key_id);
//      }
//    }
//  }
//  // Update total occupancy
//  vector_update_used(matches->rbuf_pos,out);
//  vector_update_used(matches->rbuf_mismatches,misms);
//  // Reset the intervals vector (as they had been decoded)
//  vector_clean(matches->rbuf_int);
//  return;
//}
//
//
//
//
//
//
//
//
//
//
//
//
//


