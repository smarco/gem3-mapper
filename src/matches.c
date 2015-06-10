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
  matches->total_matches_count = 0;
  matches->max_swg_score = INT32_MIN;
  // Interval Matches
  matches->interval_matches = vector_new(MATCHES_INIT_INTERVAL_MATCHES,match_interval_t);
  // Position Matches
  matches->position_matches = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t);
  matches->begin_pos_matches = ihash_new();
  matches->end_pos_matches = ihash_new();
  // CIGAR buffer
  matches->cigar_vector = vector_new(MATCHES_INIT_CIGAR_OPS,cigar_element_t);
  // Return
  return matches;
}
GEM_INLINE void matches_configure(matches_t* const matches,text_collection_t* const text_collection) {
  // Text Collection Buffer
  matches->text_collection = text_collection;
}
GEM_INLINE void matches_clear_index(matches_t* const matches) {
  ihash_clear(matches->begin_pos_matches);
  ihash_clear(matches->end_pos_matches);
}
GEM_INLINE void matches_clear(matches_t* const matches) {
  MATCHES_CHECK(matches);
  matches->max_complete_stratum = ALL;
  vector_clear(matches->counters);
  matches->min_counter_value = UINT32_MAX;
  matches->max_counter_value = 0;
  matches->total_matches_count = 0;
  matches->max_swg_score = INT32_MIN;
  vector_clear(matches->interval_matches);
  vector_clear(matches->position_matches);
  matches_clear_index(matches);
  vector_clear(matches->cigar_vector);
}
GEM_INLINE void matches_delete(matches_t* const matches) {
  MATCHES_CHECK(matches);
  // Delete all
  vector_delete(matches->counters);
  vector_delete(matches->interval_matches);
  vector_delete(matches->position_matches);
  ihash_delete(matches->begin_pos_matches);
  ihash_delete(matches->end_pos_matches);
  vector_delete(matches->cigar_vector);
  // Delete handler
  mm_free(matches);
}
/*
 * Accessors
 */
GEM_INLINE bool matches_is_mapped(const matches_t* const matches) {
  return vector_get_used(matches->position_matches) > 0 || vector_get_used(matches->interval_matches) > 0;
}
/*
 * Counters
 */
GEM_INLINE void matches_counters_calculate_min__max(matches_t* const matches) {
  // Search min/max counter
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
  // Search min/max swg-score
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  matches->max_swg_score = INT32_MIN;
  for (i=0;i<num_matches;++i) {
    if (match[i].swg_score > matches->max_swg_score) matches->max_swg_score = match[i].swg_score;
  }
}
GEM_INLINE void matches_counters_add(
    matches_t* const matches,const uint64_t distance,
    const int32_t swg_score,const uint64_t num_matches) {
  vector_t* const counters = matches->counters;
  // Reserve Memory
  if (distance >= vector_get_used(counters)) {
    vector_reserve(counters,distance+1,true);
    vector_set_used(counters,distance+1);
  }
  *vector_get_elm(counters,distance,uint64_t) += num_matches; // Add matches
  matches->total_matches_count += num_matches;
  // Update Min/Max
  matches->min_counter_value = MIN(matches->min_counter_value,distance);
  matches->max_counter_value = MAX(matches->max_counter_value,distance);
  matches->max_swg_score = MAX(matches->max_swg_score,swg_score);
}
GEM_INLINE void matches_counters_dec(
    matches_t* const matches,const uint64_t distance,const uint64_t num_matches) {
  uint64_t* const counter = vector_get_elm(matches->counters,distance,uint64_t);
  *counter -= num_matches; // Dec matches
  matches->total_matches_count -= num_matches;
  if (*counter==0) {
    matches_counters_calculate_min__max(matches); // Update Min/Max
  }
}
GEM_INLINE uint64_t matches_counters_get_min_distance(matches_t* const matches) {
  return matches->min_counter_value;
}
GEM_INLINE uint64_t matches_counters_get_max_distance(matches_t* const matches) {
  return matches->max_counter_value;
}
GEM_INLINE int32_t matches_counters_get_max_swg_score(matches_t* const matches) {
  return matches->max_swg_score;
}
GEM_INLINE uint64_t matches_counters_compact(matches_t* const matches) {
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  int64_t i = vector_get_used(matches->counters)-1;
  while (i>=0 && counters[i]==0) --i;
  vector_set_used(matches->counters,++i);
  return i;
}
GEM_INLINE uint64_t matches_counters_get_count(matches_t* const matches,const uint64_t distance) {
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  return counters[distance];
}
GEM_INLINE uint64_t matches_counters_get_total_count(matches_t* const matches) {
  return matches->total_matches_count;
}
/*
 * Matches
 */
GEM_INLINE match_trace_t* matches_get_match_traces(const matches_t* const matches) {
  return vector_get_mem(matches->position_matches,match_trace_t);
}
GEM_INLINE uint64_t matches_get_num_match_traces(const matches_t* const matches) {
  return vector_get_used(matches->position_matches);
}
GEM_INLINE void matches_get_clear_match_traces(const matches_t* const matches) {
  vector_clear(matches->position_matches);
}
GEM_INLINE cigar_element_t* match_trace_get_cigar_buffer(const matches_t* const matches,const match_trace_t* const match_trace) {
  return vector_get_elm(matches->cigar_vector,match_trace->match_alignment.cigar_offset,cigar_element_t);
}
GEM_INLINE uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace) {
  return match_trace->match_alignment.cigar_length;
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
  const cigar_element_t* cigar_element = vector_get_elm(matches->cigar_vector,cigar_buffer_offset,cigar_element_t);
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
  ihash_insert(matches->begin_pos_matches,begin_position,match_trace_offset_1);
  uint64_t* const match_trace_offset_2 = mm_stack_malloc_uint64(mm_stack);
  *match_trace_offset_2 = match_trace_offset;
  ihash_insert(matches->end_pos_matches,begin_position+effective_length,match_trace_offset_2);
}
GEM_INLINE uint64_t* matches_lookup_match(
    matches_t* const matches,const uint64_t begin_position,const uint64_t effective_length) {
  uint64_t* match_trace_offset = ihash_get(matches->begin_pos_matches,begin_position,uint64_t);
  return (match_trace_offset!=NULL) ? match_trace_offset :
      ihash_get(matches->end_pos_matches,begin_position+effective_length,uint64_t);
}
GEM_INLINE bool matches_add_match_trace_t(
    matches_t* const matches,match_trace_t* const match_trace,
    const bool update_counters,mm_stack_t* const mm_stack) {
  // Check duplicates
  uint64_t* dup_match_trace_offset = matches_lookup_match(matches,
      match_trace->match_alignment.match_position,match_trace->match_alignment.effective_length);
  if (dup_match_trace_offset==NULL) {
    // Add the match-trace
    const uint64_t new_match_trace_offset = vector_get_used(matches->position_matches);
    match_trace_t* new_match_trace;
    vector_alloc_new(matches->position_matches,match_trace_t,new_match_trace);
    *new_match_trace = *match_trace;
    // Index
    matches_index_match(matches,new_match_trace_offset,new_match_trace->match_alignment.match_position,
        new_match_trace->match_alignment.effective_length,mm_stack);
    // Update counters
    if (update_counters) matches_counters_add(matches,match_trace->distance,match_trace->swg_score,1);
    return true;
  } else {
    match_trace_t* const dup_match_trace = vector_get_elm(matches->position_matches,*dup_match_trace_offset,match_trace_t);
    if (dup_match_trace->distance > match_trace->distance) {
      // Pick up the least distant match
      if (update_counters) {
        matches_counters_dec(matches,dup_match_trace->distance,1);
        matches_counters_add(matches,match_trace->distance,match_trace->swg_score,1);
      }
      // Add the match-trace
      *dup_match_trace = *match_trace;
    }
    return false;
  }
}
GEM_INLINE void matches_add_interval_match(
    matches_t* const matches,const uint64_t lo,const uint64_t hi,
    const uint64_t length,const uint64_t distance,const bool emulated_rc_search) {
  // Check non-empty interval
  const uint64_t num_matches = hi-lo;
  if (gem_expect_false(num_matches==0)) return;
  // Setup the interval-match
  match_interval_t match_interval;
  match_interval.lo = lo;
  match_interval.hi = hi;
  match_interval.length = length;
  match_interval.distance = distance;
  match_interval.emulated_rc_search = emulated_rc_search;
  match_interval.text = NULL;
  // Update counters
  matches_counters_add(matches,distance,BOUNDED_SUBTRACTION(length,distance,0),num_matches);
  // Add the interval
  vector_insert(matches->interval_matches,match_interval,match_interval_t);
}
GEM_INLINE void matches_add_interval_set(
    matches_t* const matches,interval_set_t* const interval_set,
    const uint64_t length,const bool emulated_rc_search) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    matches_add_interval_match(matches,interval->lo,
        interval->hi,length,interval->distance,emulated_rc_search);
  }
}
GEM_INLINE void matches_hint_allocate_match_trace(matches_t* const matches,const uint64_t num_matches_trace_to_add) {
  vector_reserve_additional(matches->position_matches,num_matches_trace_to_add);
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
GEM_INLINE void matches_cigar_vector_append_insertion(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,
    const uint64_t indel_length,uint8_t* const indel_text) {
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_vector,cigar_element_t);
    if (cigar_element->type==cigar_ins) {
      cigar_element->indel.indel_length += indel_length;
      cigar_element->indel.indel_text = indel_text; // FIXME Combine indel text !!
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_ins;
  cigar_element->indel.indel_length = indel_length;
  cigar_element->indel.indel_text = indel_text;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
GEM_INLINE void matches_cigar_vector_append_deletion(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,const uint64_t indel_length) {
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_vector,cigar_element_t);
    if (cigar_element->type==cigar_del) {
      cigar_element->indel.indel_length += indel_length;
      cigar_element->indel.indel_text = NULL;
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_del;
  cigar_element->indel.indel_length = indel_length;
  cigar_element->indel.indel_text = NULL;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
GEM_INLINE void matches_cigar_vector_append_match(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,const uint64_t match_length) {
  // Check previous cigar-element (for merging)
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_vector,cigar_element_t);
    if (cigar_element->type==cigar_match) {
      cigar_element->match_length += match_length;
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_match;
  cigar_element->match_length = match_length;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
GEM_INLINE void matches_cigar_vector_append_mismatch(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,const uint8_t mismatch) {
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_mismatch;
  cigar_element->mismatch = mismatch;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
GEM_INLINE void matches_cigar_vector_append_cigar_element(
    vector_t* const cigar_vector,uint64_t* const cigar_length,cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match:
      matches_cigar_vector_append_match(cigar_vector,cigar_length,cigar_element->match_length);
      break;
    case cigar_mismatch:
      matches_cigar_vector_append_mismatch(cigar_vector,cigar_length,cigar_element->mismatch);
      break;
    case cigar_ins:
      matches_cigar_vector_append_insertion(cigar_vector,cigar_length,
          cigar_element->indel.indel_length,cigar_element->indel.indel_text);
      break;
    case cigar_del:
      matches_cigar_vector_append_deletion(cigar_vector,cigar_length,
          cigar_element->indel.indel_length);
      break;
    case cigar_soft_trim:
    case cigar_null:
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void matches_cigar_reverse(
    matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match (Even all-matching matches have CIGAR=1)
  gem_fatal_check(cigar_length==0,MATCHES_CIGAR_ZERO_LENGTH);
  // Reverse CIGAR
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_buffer_offset,cigar_element_t);
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
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_buffer_offset,cigar_element_t);
  const uint64_t middle_point = cigar_length/2;
  uint64_t i;
  for (i=0;i<middle_point;++i) {
    cigar_element_t* const origin = cigar_buffer + i;
    cigar_element_t* const flipped = cigar_buffer + (cigar_length-1-i);
    SWAP(*origin,*flipped);
  }
}
GEM_INLINE uint64_t matches_cigar_compute_edit_distance(
    const matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match
  if (cigar_length==0) return 0;
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_buffer_offset,cigar_element_t);
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
GEM_INLINE uint64_t matches_cigar_compute_edit_distance__excluding_clipping(
    const matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  // Exact Match
  if (cigar_length==0) return 0;
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_buffer_offset,cigar_element_t);
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
GEM_INLINE int64_t matches_cigar_element_effective_length(const cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match:
      return cigar_element->match_length;
      break;
    case cigar_mismatch:
      return 1;
      break;
    case cigar_del:
    case cigar_soft_trim:
      break;
    case cigar_ins:
      return cigar_element->indel.indel_length;
      break;
    case cigar_null:
    default:
      GEM_INVALID_CASE();
      break;
  }
  return 0;
}
GEM_INLINE int64_t matches_cigar_effective_length(
    vector_t* const cigar_vector,const uint64_t cigar_offset,const uint64_t cigar_length) {
  uint64_t effective_length = 0;
  // Traverse all CIGAR elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t i;
  for (i=0;i<cigar_length;++i) {
    effective_length += matches_cigar_element_effective_length(cigar_buffer+i);
  }
  // Return effective length
  return effective_length;
}
/*
 * Sorting Matches
 */
int match_trace_cmp_distance(const match_trace_t* const a,const match_trace_t* const b) {
  return (int)a->distance - (int)b->distance;
}
GEM_INLINE void matches_sort_by_distance(matches_t* const matches) {
  qsort(vector_get_mem(matches->position_matches,match_trace_t),
      vector_get_used(matches->position_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_distance);
}
int match_trace_cmp_swg_score(const match_trace_t* const a,const match_trace_t* const b) {
  return (int)b->swg_score - (int)a->swg_score;
}
GEM_INLINE void matches_sort_by_swg_score(matches_t* const matches) {
  qsort(vector_get_mem(matches->position_matches,match_trace_t),
      vector_get_used(matches->position_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_swg_score);
}
int match_trace_cmp_mapq_score(const match_trace_t* const a,const match_trace_t* const b) {
  return (int)b->mapq_score - (int)a->mapq_score;
}
GEM_INLINE void matches_sort_by_mapq_score(matches_t* const matches) {
  qsort(vector_get_mem(matches->position_matches,match_trace_t),
      vector_get_used(matches->position_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_mapq_score);
}
int match_trace_cmp_sequence_name__position(const match_trace_t* const a,const match_trace_t* const b) {
  const int cmp = gem_strcmp(a->sequence_name,b->sequence_name);
  return (cmp!=0) ? cmp : ((int)a->text_position - (int)b->text_position);
}
GEM_INLINE void matches_sort_by_sequence_name__position(matches_t* const matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(matches->position_matches,match_trace_t),
      vector_get_used(matches->position_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_sequence_name__position);
}
/*
 * Curation
 */
GEM_INLINE void matches_curate(matches_t* const matches,const double swg_score_difference) {
  uint64_t* const counter = vector_get_mem(matches->counters,uint64_t);
  match_trace_t* match_in = matches_get_match_traces(matches);
  match_trace_t* match_out = match_in;
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Delete matches with bigger error than the best
  const int32_t min_swg_score = (swg_score_difference >= 1.0) ? (int32_t) swg_score_difference :
      matches->max_swg_score - (int32_t)(swg_score_difference*(double)matches->max_swg_score);
  uint64_t i;
  for (i=0;i<num_matches;++i,++match_in) {
    if (match_in->swg_score < min_swg_score) {
      --(counter[match_in->distance]);
      --(matches->total_matches_count);
    } else {
      if (match_out != match_in) *match_out = *match_in;
      ++match_out;
    }
  }
  vector_update_used(matches->position_matches,match_out);
  // Recompact counters
  matches_counters_compact(matches);
  // Recalculate min/max
  matches_counters_calculate_min__max(matches);
}
/*
 * Logistic Regression
 */
GEM_INLINE void matches_metrics_compute(
    matches_t* const matches,const uint64_t read_length,
    double* const id,double* const sub_id,
    double* const swg_ratio,double* const sub_swg_ratio,
    uint64_t* const mcs,uint64_t* const fs_matches,uint64_t* const sub_matches) {
  // Compute identities
  const uint64_t num_counters = vector_get_used(matches->counters);
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  uint64_t i=0, j;
  while (i<num_counters && counters[i]==0) ++i;
  if (i == num_counters) {
    *id = 0.0;
    *sub_id = 0.0;
  } else {
    *id = (double)(read_length-i)/(double)read_length;
    j=i+1;
    while (j<num_counters && counters[j]==0) ++j;
    if (j == num_counters) {
      *sub_id = 0.0;
    } else {
      *sub_id = (double)(read_length-j)/(double)read_length;
    }
  }
  // Compute swg
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  int32_t max1_swg = 0, max2_swg = 0;
  for (i=0;i<num_matches;++i) {
    if (match[i].swg_score > max1_swg) {
      if (max1_swg > max2_swg) max2_swg = max1_swg;
      max1_swg = match[i].swg_score;
    } else if (match[i].swg_score > max2_swg) {
      max2_swg = match[i].swg_score;
    }
  }
  *swg_ratio = (double)max1_swg/(double)read_length;
  *sub_swg_ratio = (double)max2_swg/(double)read_length;
  // Others
  *mcs = matches->max_complete_stratum;
  *fs_matches = matches_counters_get_count(matches,matches_counters_get_min_distance(matches));
  *sub_matches = matches_counters_get_total_count(matches) - *fs_matches;
}
GEM_INLINE double matches_classify_unique(matches_t* const matches,const uint64_t mcs) {
  /*
   * Trivially discards mmaps & ties and returns the probability of
   * the first position-match (primary match) of being a true positive
   */
  // No matches
  if (matches_counters_get_total_count(matches)==0) return 0.0;
  // Metrics
  double id, sub_id, swg_ratio, sub_swg_ratio;
  uint64_t fs_matches, sub_matches, dummy;
  matches_metrics_compute(matches,100,
      &id,&sub_id,&swg_ratio,&sub_swg_ratio,&dummy,&fs_matches,&sub_matches);
  // Remove ties & mmaps
  if (fs_matches > 1 || sub_matches > 0) return 0.0;
  // Compute PR
  const double lr_factor = -31.886 + swg_ratio * 29.617 + (double)mcs * 4.695;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
GEM_INLINE double matches_classify_ambiguous(matches_t* const matches,const uint64_t mcs) {
  /*
   * Trivially discards ties and returns the (1-probability) of
   * the first position-match (primary match) of being a false positive (ambiguous map)
   */
  // No matches
  if (matches_counters_get_total_count(matches)==0) return 0.0;
  // Metrics
  double id, sub_id, swg_ratio, sub_swg_ratio;
  uint64_t fs_matches, sub_matches, dummy;
  matches_metrics_compute(matches,100,
      &id,&sub_id,&swg_ratio,&sub_swg_ratio,&dummy,&fs_matches,&sub_matches);
  // Remove ties
  if (fs_matches > 1) return 0.0;
  // Compute PR
  const double lr_factor = -18.73331 + sub_id * 0.81486 + swg_ratio * 19.93329 +
      sub_swg_ratio * -2.73863 + (double)mcs * 2.30287 + (double)sub_matches * -0.01749;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
GEM_INLINE double matches_classify_mmaps(matches_t* const matches,const uint64_t mcs) {
  /*
   * Trivially discards ties and returns the probability of
   * the first position-match (primary match) of being a true positive
   */
  // No matches
  if (matches_counters_get_total_count(matches)==0) return 0.0;
  // Metrics
  double id, sub_id, swg_ratio, sub_swg_ratio;
  uint64_t fs_matches, sub_matches, dummy;
  matches_metrics_compute(matches,100,
      &id,&sub_id,&swg_ratio,&sub_swg_ratio,&dummy,&fs_matches,&sub_matches);
  // Remove ties
  if (fs_matches > 1) return 0.0;
  // Mmaps
  const double lr_factor = 16.261 +
      (double)id * -134.234 +
      (double)sub_id * 102.120 +
      (double)swg_ratio * 77.025 +
      (double)sub_swg_ratio * -59.739 +
      (double)mcs * 0.452;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
GEM_INLINE double matches_classify_ties(matches_t* const matches) {
  /*
   * Trivially discards mmaps & ties and returns the probability of
   * the first position-match (primary match) of being a true positive
   */
  // No matches
  if (matches_counters_get_total_count(matches)==0) return 0.0;
  const uint64_t fs_matches = matches_counters_get_count(matches,matches_counters_get_min_distance(matches));
  if (fs_matches > 1) return 1.0;
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches > 1 && match[0].swg_score==match[1].swg_score) return 1.0;
  return 0.0;
}
//GEM_INLINE bool matches_resolved_as_signal(matches_t* const matches,const uint64_t mcs) {
//  // Isolate Pure-Signal
//  const double pr = matches_logit_psignal(matches,mcs);
//  if (pr >= 0.98) return true;
//  // Return unknown
//  return false;
//}
//GEM_INLINE bool matches_resolved_as_noise(matches_t* const matches,const uint64_t mcs) {
//  // Isolate Pure-Noise
//  const double pr = matches_logit_pnoise(matches,mcs);
//  if (pr <= 0.20) return true;
//  // Return unknown
//  return false;
//}
/*
 * Display
 */
GEM_INLINE void match_cigar_print(
    FILE* const stream,vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  uint64_t j;
  for (j=0;j<cigar_length;++j) {
    cigar_element_t* const cigar_element = vector_get_elm(cigar_vector,cigar_buffer_offset+j,cigar_element_t);
    switch (cigar_element->type) {
      case cigar_match: tab_fprintf(stream,"%luM",cigar_element->match_length); break;
      case cigar_mismatch: tab_fprintf(stream,"1X"); break;
      case cigar_ins: tab_fprintf(stream,"%luI",cigar_element->indel.indel_length); break;
      case cigar_del: tab_fprintf(stream,"%luD",cigar_element->indel.indel_length); break;
      default: GEM_INVALID_CASE(); break;
    }
  }
}
GEM_INLINE void matches_metrics_print(matches_t* const matches) {
  // Metrics
  double id, sub_id, swg_ratio, sub_swg_ratio;
  uint64_t mcs, fs_matches, sub_matches;
  matches_metrics_compute(matches,100,
      &id,&sub_id,&swg_ratio,&sub_swg_ratio,&mcs,&fs_matches,&sub_matches);
  // id
  fprintf(stdout,"%f\t",id);
  // sub_id
  fprintf(stdout,"%f\t",sub_id);
  // swg
  fprintf(stdout,"%f\t",swg_ratio);
  // sub_swg
  fprintf(stdout,"%f\t",sub_swg_ratio);
  // swg_diff
  fprintf(stdout,"%f\t",swg_ratio-sub_swg_ratio);
  // mcs
  fprintf(stdout,"%lu\t",mcs);
  // fs_matches
  fprintf(stdout,"%lu\t",fs_matches);
  // sub_matches
  fprintf(stdout,"%lu\t",sub_matches);
  // mapq_score
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches > 0) {
    match_trace_t* const match = matches_get_match_traces(matches);
    fprintf(stdout,"%d\n",match[0].mapq_score);
  } else {
    fprintf(stdout,"\n");
  }
}

