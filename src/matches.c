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
  // Search-matches state
  matches->max_complete_stratum = ALL;
  // Matches Counters
  matches->counters = vector_new(MATCHES_INIT_COUNTERS,uint64_t);
  matches->min_counter_value = UINT32_MAX;
  matches->max_counter_value = 0;
  // Interval Matches
  matches->interval_matches = vector_new(MATCHES_INIT_INTERVAL_MATCHES,match_interval_t);
  // Position Matches // TODO IF vectors -> iterator for HUGE datasets
  matches->global_matches = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t);
  matches->begin_gmatches = ihash_new();
  matches->end_gmatches = ihash_new();
  // CIGAR buffer
  matches->cigar_buffer = vector_new(MATCHES_INIT_CIGAR_OPS,cigar_element_t);
  // Restore Point (RP)
  matches->rp_counters = vector_new(MATCHES_INIT_COUNTERS,uint64_t);
  matches->rp_interval_matches_used = 0;
  matches->rp_global_matches_used = 0;
  matches->rp_cigar_buffer_used = 0;
  // Return
  return matches;
}
GEM_INLINE void matches_configure(matches_t* const matches,text_collection_t* const text_collection) {
  // Text Collection Buffer
  matches->text_collection = text_collection;
}
GEM_INLINE void matches_clear(matches_t* const matches) {
  MATCHES_CHECK(matches);
  matches->max_complete_stratum = ALL;
  vector_clear(matches->counters);
  matches->min_counter_value = UINT32_MAX;
  matches->max_counter_value = 0;
  vector_clear(matches->interval_matches);
  vector_clear(matches->global_matches);
  ihash_clear(matches->begin_gmatches);
  ihash_clear(matches->end_gmatches);
  vector_clear(matches->cigar_buffer);
}
GEM_INLINE void matches_delete(matches_t* const matches) {
  MATCHES_CHECK(matches);
  // Delete all
  vector_delete(matches->counters);
  vector_delete(matches->interval_matches);
  vector_delete(matches->global_matches);
  ihash_delete(matches->begin_gmatches);
  ihash_delete(matches->end_gmatches);
  vector_delete(matches->cigar_buffer);
  vector_delete(matches->rp_counters);
  // Delete handler
  mm_free(matches);
}
/*
 * Counters
 */
GEM_INLINE void matches_counters_calculate_min__max(matches_t* const matches) {
  const uint64_t num_counters = vector_get_used(matches->counters);
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  uint64_t i = 0;
  while (i<num_counters && counters[i]==0) ++i;
  if (i==num_counters) {
    matches->min_counter_value = UINT32_MAX;
    matches->max_counter_value = 0;
  } else {
    matches->min_counter_value = i;
    uint64_t max = i;
    while (i<num_counters) {
      if (counters[i]!=0) max = i;
      ++i;
    }
    matches->max_counter_value = max;
  }
}
GEM_INLINE void matches_counters_add(matches_t* const matches,const uint64_t distance,const uint64_t num_matches) {
  vector_t* const counters = matches->counters;
  // Reserve Memory
  if (distance >= vector_get_used(counters)) {
    vector_reserve(counters,distance+1,true);
    vector_set_used(counters,distance+1);
  }
  *vector_get_elm(counters,distance,uint64_t) += num_matches; // Add matches
  // Update Min/Max
  matches->min_counter_value = MIN(matches->min_counter_value,distance);
  matches->max_counter_value = MAX(matches->max_counter_value,distance);
}
GEM_INLINE void matches_counters_dec(matches_t* const matches,const uint64_t distance,const uint64_t num_matches) {
  uint64_t* counter = vector_get_elm(matches->counters,distance,uint64_t);
  *counter -= num_matches; // Dec matches
  if (*counter==0) {
    matches_counters_calculate_min__max(matches); // Update Min/Max
  }
}
GEM_INLINE uint64_t matches_counters_get_min(matches_t* const matches) {
  return matches->min_counter_value;
}
GEM_INLINE uint64_t matches_counters_get_max(matches_t* const matches) {
  return matches->max_counter_value;
}
GEM_INLINE uint64_t matches_counters_compact(matches_t* const matches) {
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  int64_t i = vector_get_used(matches->counters)-1;
  while (i>=0 && counters[i]==0) --i;
  vector_set_used(matches->counters,++i);
  return i;
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
GEM_INLINE uint64_t match_trace_get_distance(const match_trace_t* const match_trace) {
  return match_trace->distance;
}
GEM_INLINE int64_t match_trace_get_effective_length(
    matches_t* const matches,const uint64_t read_length,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match
  if (cigar_length==0) return read_length; // Even all-matching matches have CIGAR=1
  // Traverse CIGAR
  const cigar_element_t* cigar_element = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  int64_t i, effective_length = read_length;
  for (i=0;i<cigar_length;++i,++cigar_element) {
    switch (cigar_element->type) {
      case cigar_ins:
        effective_length += cigar_element->indel.indel_length;
        break;
      case cigar_del:
      case cigar_soft_trim:
        effective_length -= cigar_element->indel.indel_length;
        break;
      default:
        break;
    }
  }
  GEM_INTERNAL_CHECK(effective_length >= 0,"Match effective length must be positive");
  return effective_length;
}
/*
 * Adding Matches
 */
GEM_INLINE void matches_index_match(
    matches_t* const matches, const uint64_t match_trace_offset,
    const uint64_t begin_position,const uint64_t effective_length,
    mm_stack_t* const mm_stack) {
  // Store begin/end positions of the match (as to fast index matches)
  uint64_t* const match_trace_offset_1 = mm_stack_malloc_uint64(mm_stack);
  *match_trace_offset_1 = match_trace_offset;
  ihash_insert(matches->begin_gmatches,begin_position,match_trace_offset_1);
  uint64_t* const match_trace_offset_2 = mm_stack_malloc_uint64(mm_stack);
  *match_trace_offset_2 = match_trace_offset;
  ihash_insert(matches->end_gmatches,begin_position+effective_length,match_trace_offset_2);
}
GEM_INLINE uint64_t* matches_lookup_match(
    matches_t* const matches,const uint64_t begin_position,const uint64_t effective_length) {
  uint64_t* match_trace_offset = ihash_get(matches->begin_gmatches,begin_position,uint64_t);
  return (match_trace_offset!=NULL) ? match_trace_offset :
      ihash_get(matches->end_gmatches,begin_position+effective_length,uint64_t);
}
GEM_INLINE void matches_add_match_trace_t(
    matches_t* const matches,match_trace_t* const match_trace,
    const bool update_counters,mm_stack_t* const mm_stack) {
  // Check duplicates
  uint64_t* dup_match_trace_offset = matches_lookup_match(matches,match_trace->index_position,match_trace->effective_length);
  if (dup_match_trace_offset==NULL) {
    // Add the match-trace
    const uint64_t new_match_trace_offset = vector_get_used(matches->global_matches);
    match_trace_t* new_match_trace;
    vector_alloc_new(matches->global_matches,match_trace_t,new_match_trace);
    *new_match_trace = *match_trace;
    // Index
    matches_index_match(matches,new_match_trace_offset,
        new_match_trace->index_position,new_match_trace->effective_length,mm_stack);
    // Update counters
    if (update_counters) matches_counters_add(matches,match_trace->distance,1);
  } else {
    match_trace_t* const dup_match_trace = vector_get_elm(matches->global_matches,*dup_match_trace_offset,match_trace_t);
    if (dup_match_trace->distance > match_trace->distance) {
      // Pick up the least distant match
      if (update_counters) {
        matches_counters_dec(matches,dup_match_trace->distance,1);
        matches_counters_add(matches,match_trace->distance,1);
      }
      // Add the match-trace
      *dup_match_trace = *match_trace;
    }
  }
}
GEM_INLINE void matches_add_interval_match(
    matches_t* const matches,const uint64_t lo,const uint64_t hi,
    const uint64_t length,const uint64_t distance,const strand_t strand) {
  // Check non-empty interval
  const uint64_t num_matches = hi-lo;
  if (gem_expect_false(num_matches==0)) return;
  // Setup the interval-match
  match_interval_t match_interval;
  match_interval.lo = lo;
  match_interval.hi = hi;
  match_interval.length = length;
  match_interval.distance = distance;
  match_interval.strand = strand;
  match_interval.text = NULL;
  // Update counters
  matches_counters_add(matches,distance,num_matches);
  // Add the interval
  vector_insert(matches->interval_matches,match_interval,match_interval_t);
}
GEM_INLINE void matches_add_interval_set(
    matches_t* const matches,interval_set_t* const interval_set,
    const uint64_t length,const strand_t strand) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    matches_add_interval_match(matches,interval->lo,
        interval->hi,length,interval->distance,strand);
  }
}
GEM_INLINE void matches_hint_allocate_match_trace(matches_t* const matches,const uint64_t num_matches_trace_to_add) {
  vector_reserve_additional(matches->global_matches,num_matches_trace_to_add);
}
GEM_INLINE void matches_hint_allocate_match_interval(matches_t* const matches,const uint64_t num_matches_interval_to_add) {
  vector_reserve_additional(matches->interval_matches,num_matches_interval_to_add);
}
/*
 * CIGAR Handling
 */
GEM_INLINE void matches_cigar_buffer_add_cigar_element(
    cigar_element_t** const cigar_buffer_sentinel,const cigar_t cigar_element_type,
    const uint64_t element_length,uint8_t* const indel_text) {
  if ((*cigar_buffer_sentinel)->type == cigar_element_type) {
    if ((*cigar_buffer_sentinel)->type == cigar_match) {
      (*cigar_buffer_sentinel)->match_length += element_length;
    } else {
      (*cigar_buffer_sentinel)->indel.indel_length += element_length;
      if ((*cigar_buffer_sentinel)->type == cigar_ins) {
        (*cigar_buffer_sentinel)->indel.indel_text = indel_text;
      }
    }
  } else {
    if ((*cigar_buffer_sentinel)->type != cigar_null) ++(*cigar_buffer_sentinel);
    (*cigar_buffer_sentinel)->type = cigar_element_type;
    if ((*cigar_buffer_sentinel)->type == cigar_match) {
      (*cigar_buffer_sentinel)->match_length = element_length;
    } else {
      (*cigar_buffer_sentinel)->indel.indel_length = element_length;
      if ((*cigar_buffer_sentinel)->type == cigar_ins) {
        (*cigar_buffer_sentinel)->indel.indel_text = indel_text;
      }
    }
  }
}
GEM_INLINE void matches_cigar_buffer_append_indel(
    vector_t* const cigar_buffer,uint64_t* const current_cigar_length,
    const cigar_t cigar_element_type,const uint64_t element_length,uint8_t* const indel_text) {
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_buffer,cigar_element_t);
    if (cigar_element->type==cigar_element_type) {
      cigar_element->indel.indel_length += element_length;
      cigar_element->indel.indel_text = indel_text;
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_buffer,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_buffer,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_element_type;
  cigar_element->indel.indel_length = element_length;
  cigar_element->indel.indel_text = indel_text;
  vector_inc_used(cigar_buffer); // Increment used
  *current_cigar_length += 1;
}
GEM_INLINE void matches_cigar_buffer_append_match(
    vector_t* const cigar_buffer,uint64_t* const current_cigar_length,
    const uint64_t match_length) {
  // Check previous cigar-element (for merging)
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_buffer,cigar_element_t);
    if (cigar_element->type==cigar_match) {
      cigar_element->match_length += match_length;
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_buffer,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_buffer,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_match;
  cigar_element->match_length = match_length;
  vector_inc_used(cigar_buffer); // Increment used
  *current_cigar_length += 1;
}
GEM_INLINE void matches_cigar_buffer_append_mismatch(
    vector_t* const cigar_buffer,uint64_t* const current_cigar_length,
    const cigar_t cigar_element_type,const uint8_t mismatch) {
  vector_reserve_additional(cigar_buffer,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_buffer,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_element_type;
  cigar_element->mismatch = mismatch;
  vector_inc_used(cigar_buffer); // Increment used
  *current_cigar_length += 1;
}
GEM_INLINE void matches_cigar_reverse(
    matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match (Even all-matching matches have CIGAR=1)
  gem_fatal_check(cigar_length==0,MATCHES_CIGAR_ZERO_LENGTH);
  // Reverse CIGAR
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  const uint64_t middle_point = cigar_length/2;
  uint64_t i;
  for (i=0;i<middle_point;++i) {
    cigar_element_t* const origin = cigar_buffer + i;
    cigar_element_t* const flipped = cigar_buffer + (cigar_length-1-i);
    SWAP(*origin,*flipped);
    if (origin->type == cigar_mismatch) origin->mismatch = dna_encoded_complement(origin->mismatch);
    if (flipped->type == cigar_mismatch) flipped->mismatch = dna_encoded_complement(flipped->mismatch);
    // FIXME In case of indel, flip @origin->indel.indel_text (only SAM.MD field is using it)
  }
  if (cigar_length%2) {
    cigar_element_t* const middle = cigar_buffer + middle_point;
    if (middle->type == cigar_mismatch) middle->mismatch = dna_encoded_complement(middle->mismatch);
    // FIXME In case of indel, flip @origin->indel.indel_text (only SAM.MD field is using it)
  }
}
GEM_INLINE void matches_cigar_reverse_colorspace(
    matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match
  if (cigar_length==0) return; // Even all-matching matches have CIGAR=1
  // Reverse CIGAR
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  const uint64_t middle_point = cigar_length/2;
  uint64_t i;
  for (i=0;i<middle_point;++i) {
    cigar_element_t* const origin = cigar_buffer + i;
    cigar_element_t* const flipped = cigar_buffer + (cigar_length-1-i);
    SWAP(*origin,*flipped);
  }
}
GEM_INLINE uint64_t matches_cigar_calculate_edit_distance(
    const matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match
  if (cigar_length==0) return 0;
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  uint64_t i, edit_distance = 0;
  for (i=0;i<cigar_length;++i) {
    switch (cigar_buffer[i].type) {
      case cigar_match: break;
      case cigar_mismatch:
        edit_distance += cigar_buffer[i].type;
        break;
      case cigar_ins:
      case cigar_del:
      case cigar_soft_trim:
        edit_distance += cigar_buffer[i].indel.indel_length;
        break;
      case cigar_null:
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  return edit_distance;
}
GEM_INLINE uint64_t matches_cigar_calculate_edit_distance__excluding_clipping(
    const matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match
  if (cigar_length==0) return 0;
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_buffer,cigar_buffer_offset,cigar_element_t);
  uint64_t i, edit_distance = 0;
  for (i=0;i<cigar_length;++i) {
    switch (cigar_buffer[i].type) {
      case cigar_match: break;
      case cigar_mismatch:
        edit_distance += cigar_buffer[i].type;
        break;
      case cigar_del:
      case cigar_soft_trim:
        if (i==0 || i==cigar_length-1) break;
      // No break
      case cigar_ins:
        edit_distance += cigar_buffer[i].indel.indel_length;
        break;
      case cigar_null:
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  return edit_distance;
}
/*
 * Status
 */
GEM_INLINE bool matches_is_mapped(const matches_t* const matches) {
  return vector_get_used(matches->global_matches) > 0 || vector_get_used(matches->interval_matches) > 0;
}
/*
 * Sorting Matches
 */
int match_trace_cmp_distance(const match_trace_t* const a,const match_trace_t* const b) {
  return (int)a->distance - (int)b->distance;
}
GEM_INLINE void matches_sort_by_distance(matches_t* const matches) {
  qsort(vector_get_mem(matches->global_matches,match_trace_t),
      vector_get_used(matches->global_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_distance);
}
int match_trace_cmp_mapq_score(const match_trace_t* const a,const match_trace_t* const b) {
  return (int)b->mapq_score - (int)a->mapq_score;
}
GEM_INLINE void matches_sort_by_mapq_score(matches_t* const matches) {
  qsort(vector_get_mem(matches->global_matches,match_trace_t),
      vector_get_used(matches->global_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_mapq_score);
}
int match_trace_cmp_sequence_name__position(const match_trace_t* const a,const match_trace_t* const b) {
  const int cmp = gem_strcmp(a->sequence_name,b->sequence_name);
  return (cmp!=0) ? cmp : ((int)a->text_position - (int)b->text_position);
}
GEM_INLINE void matches_sort_by_sequence_name__position(matches_t* const matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(matches->global_matches,match_trace_t),
      vector_get_used(matches->global_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_sequence_name__position);
}


