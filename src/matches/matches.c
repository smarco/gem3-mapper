/*
 * PROJECT: GEMMapper
 * FILE: matches.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure to store alignment matches {sequence,position,strand,CIGAR}
 */

#include "matches/matches.h"
#include "matches/matches_cigar.h"
#include "matches/matches_classify.h"

/*
 * Constants
 */
#define MATCHES_INIT_INTERVAL_MATCHES  100
#define MATCHES_INIT_GLOBAL_MATCHES   1000
#define MATCHES_INIT_CIGAR_OPS        5000

/*
 * Matches Classes
 */
const char* matches_class_label[] =
{
    [0] = "unmapped",
    [1] = "tie-d0",
    [2] = "tie-d1",
    [3] = "mmap",
    [4] = "unique"
};

/*
 * Setup
 */
matches_t* matches_new() {
  // Allocate handler
  matches_t* const matches = mm_alloc(matches_t);
  // Search-matches state
  matches->matches_class = matches_class_unmapped;
  matches->max_complete_stratum = ALL;
  // Matches Counters
  matches->counters = matches_counters_new();
  // Position Matches
  matches->position_matches = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t);
  matches->begin_pos_matches = ihash_new();
  matches->end_pos_matches = ihash_new();
  // Local Matches
  matches->local_matches = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t);
  // CIGAR buffer
  matches->cigar_vector = vector_new(MATCHES_INIT_CIGAR_OPS,cigar_element_t);
  // Init metrics
  matches_metrics_init(&matches->metrics);
  // Return
  return matches;
}
void matches_configure(matches_t* const matches,text_collection_t* const text_collection) {
  // Text Collection Buffer
  matches->text_collection = text_collection;
}
void matches_clear(matches_t* const matches) {
  matches->matches_class = matches_class_unmapped;
  matches->max_complete_stratum = ALL;
  matches_counters_clear(matches->counters);
  matches_metrics_init(&matches->metrics);
  vector_clear(matches->position_matches);
  matches_index_clear(matches);
  vector_clear(matches->local_matches);
  vector_clear(matches->cigar_vector);
}
void matches_delete(matches_t* const matches) {
  // Delete all
  matches_counters_delete(matches->counters);
  vector_delete(matches->position_matches);
  ihash_delete(matches->begin_pos_matches);
  ihash_delete(matches->end_pos_matches);
  vector_delete(matches->local_matches);
  vector_delete(matches->cigar_vector);
  // Delete handler
  mm_free(matches);
}
/*
 * Accessors
 */
bool matches_is_mapped(const matches_t* const matches) {
  return vector_get_used(matches->position_matches) > 0;
}
void matches_recompute_metrics(matches_t* const matches) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  // ReCompute distance metrics
  const match_trace_t* match = matches_get_match_trace_buffer(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  matches_metrics_init(metrics);
  uint64_t i;
  for (i=0;i<num_matches;++i,++match) {
    matches_metrics_update(metrics,match->distance,match->edit_distance,match->swg_score);
  }
}
uint64_t matches_get_first_stratum_matches(matches_t* const matches) {
  const uint64_t min_distance = matches_metrics_get_min_distance(&matches->metrics);
  return (min_distance==UINT32_MAX) ? 0 : matches_counters_get_count(matches->counters,min_distance);
}
uint64_t matches_get_subdominant_stratum_matches(matches_t* const matches) {
  const uint64_t first_stratum_matches = matches_get_first_stratum_matches(matches);
  return matches_counters_get_total_count(matches->counters) - first_stratum_matches;
}
uint8_t matches_get_primary_mapq(matches_t* const matches) {
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches == 0) return 0;
  const match_trace_t* match = matches_get_match_trace_buffer(matches);
  return match->mapq_score;
}
/*
 * Index
 */
void matches_index_match(
    matches_t* const matches,
    match_trace_t* const match_trace,
    const uint64_t match_trace_offset,
    mm_stack_t* const mm_stack) {
  // Compute dimensions
  const uint64_t begin_position = match_trace->match_alignment.match_position;
  const uint64_t effective_length = match_trace->match_alignment.effective_length;
  // Allocate pointer for offset
  match_trace->match_trace_offset = mm_stack_malloc_uint64(mm_stack);
  *(match_trace->match_trace_offset) = match_trace_offset;
  // Store begin/end positions of the match (as to fast index matches)
  ihash_insert(matches->begin_pos_matches,begin_position,match_trace->match_trace_offset);
  ihash_insert(matches->end_pos_matches,begin_position+effective_length,match_trace->match_trace_offset);
}
uint64_t* matches_lookup_match(
    matches_t* const matches,
    const uint64_t begin_position,
    const uint64_t effective_length) {
  uint64_t* match_trace_offset = ihash_get(matches->begin_pos_matches,begin_position,uint64_t);
  return (match_trace_offset!=NULL) ? match_trace_offset :
      ihash_get(matches->end_pos_matches,begin_position+effective_length,uint64_t);
}
void matches_index_clear(matches_t* const matches) {
  ihash_clear(matches->begin_pos_matches);
  ihash_clear(matches->end_pos_matches);
}
void matches_index_rebuild(matches_t* const matches,mm_stack_t* const mm_stack) {
  // Clear
  matches_index_clear(matches);
  // Re-index positions
  match_trace_t* const match = matches_get_match_trace_buffer(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  uint64_t i;
  for (i=0;i<num_matches;++i) {
    matches_index_match(matches,match+i,i,mm_stack);
  }
}
/*
 * Matches
 */
match_trace_t* matches_get_match_trace_buffer(const matches_t* const matches) {
  return vector_get_mem(matches->position_matches,match_trace_t);
}
match_trace_t* matches_get_match_trace(const matches_t* const matches,const uint64_t offset) {
  return vector_get_elm(matches->position_matches,offset,match_trace_t);
}
uint64_t matches_get_num_match_traces(const matches_t* const matches) {
  return vector_get_used(matches->position_matches);
}
void matches_get_clear_match_traces(const matches_t* const matches) {
  vector_clear(matches->position_matches);
}
cigar_element_t* match_trace_get_cigar_buffer(const matches_t* const matches,const match_trace_t* const match_trace) {
  return vector_get_elm(matches->cigar_vector,match_trace->match_alignment.cigar_offset,cigar_element_t);
}
uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace) {
  return match_trace->match_alignment.cigar_length;
}
uint64_t match_trace_get_event_distance(const match_trace_t* const match_trace) {
  return match_trace->distance;
}
int64_t match_trace_get_effective_length(
    matches_t* const matches,
    const uint64_t read_length,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  // Exact Match
  if (cigar_length==0) return read_length; // Even all-matching matches have CIGAR=1
  // Traverse CIGAR
  const cigar_element_t* cigar_element = vector_get_elm(matches->cigar_vector,cigar_buffer_offset,cigar_element_t);
  int64_t i, effective_length = read_length;
  for (i=0;i<cigar_length;++i,++cigar_element) {
    switch (cigar_element->type) {
      case cigar_ins:
        effective_length += cigar_element->length;
        break;
      case cigar_del:
        effective_length -= cigar_element->length;
        break;
      default:
        break;
    }
  }
  GEM_INTERNAL_CHECK(effective_length >= 0,"Match effective length must be positive");
  return effective_length;
}
/*
 * Matches Rank Consistency
 *   Preserves the first @max_reported_matches with the least distance
 */
match_trace_t* matches_sort_preserve_rank_consistency(
    matches_t* const matches,
    select_parameters_t* const select_parameters,
    const alignment_model_t alignment_model,
    const uint64_t new_match_trace_offset) {
  // Parameters
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  match_trace_t* const match_trace = matches_get_match_trace_buffer(matches);
  match_trace_t* const new_match_trace = match_trace + new_match_trace_offset;
  if (num_matches==0 || alignment_model==alignment_model_none) {
    return new_match_trace; // No swap, return current position
  }
  // Find proper position (as to swap with newest one & preserve rank)
  uint64_t match_offset = 0;
  uint64_t current_stratum = 0, current_stratum_distance = match_trace->distance;
  while (match_offset < num_matches) {
    // Check current match distance
    const bool is_within_range = (alignment_model == alignment_model_gap_affine) ?
        match_trace[match_offset].swg_score < new_match_trace->swg_score :
        match_trace[match_offset].distance > new_match_trace->distance;
    if (is_within_range) {
      // New match-trace is within the eligible range; swap to maintain rank consistency
      match_trace_t* const current_match_trace = match_trace + match_offset;
      if (*(new_match_trace->match_trace_offset) != match_offset) {
        SWAP(*current_match_trace,*new_match_trace);
        SWAP(*(current_match_trace->match_trace_offset),*(new_match_trace->match_trace_offset));
      }
      return current_match_trace; // Return new position
    }
    // Check stratum
    if (match_trace[match_offset].distance > current_stratum_distance) {
      current_stratum_distance = match_trace[match_offset].distance;
      ++current_stratum;
    }
    // Check select limits
    if (current_stratum >= select_parameters->min_reported_strata_nominal &&
        match_offset+1 >= select_parameters->max_reported_matches) {
      // New match-trace is out of the eligible range
      return new_match_trace; // No swap, return current position
    }
    // Next
    ++match_offset;
  }
  // No swap, return current position
  return new_match_trace;
}
match_trace_t* matches_get_ranked_match_trace(
    matches_t* const matches,
    select_parameters_t* const select_parameters) {
  // Parameters
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  match_trace_t* const match_trace = matches_get_match_trace_buffer(matches);
  // Return selected match-trace from ranked ranged
  if (select_parameters->min_reported_strata_nominal == 0) {
    return match_trace + (MIN(select_parameters->max_reported_matches,num_matches)-1);
  } else {
    uint64_t match_offset = 0;
    uint64_t current_stratum = 0, current_stratum_distance = match_trace->distance;
    while (match_offset < num_matches) {
      // Check stratum
      if (match_trace[match_offset].distance > current_stratum_distance) {
        current_stratum_distance = match_trace[match_offset].distance;
        ++current_stratum;
      }
      // Check select limits
      if (current_stratum >= select_parameters->min_reported_strata_nominal &&
          match_offset+1 >= select_parameters->max_reported_matches) {
        return match_trace + match_offset; // Return last match-trace in the eligible range
      }
      // Next
      ++match_offset;
    }
    return match_trace + (num_matches-1);
  }
}
/*
 * Add Matches
 */
void match_trace_locate(match_trace_t* const match_trace,const locator_t* const locator) {
  GEM_INTERNAL_CHECK(match_trace->match_alignment.effective_length >= 0,"Match effective length must be positive");
  location_t location;
  locator_map(locator,match_trace->match_alignment.match_position,&location);
  match_trace->text_position = location.position;
  match_trace->sequence_name = location.tag;
  if (location.strand == Reverse) { // Adjust position by the effective length
    match_trace->text_position -= (uint64_t) match_trace->match_alignment.effective_length;
    match_trace->strand = Reverse;
    GEM_INTERNAL_CHECK(!match_trace->emulated_rc_search,
        "Archive-Select. Locating match-trace. "
        "Impossible combination (search_strand==Reverse while emulated searching FR-index)");
  } else {
    match_trace->strand = match_trace->emulated_rc_search ? Reverse : Forward;
  }
  match_trace->bs_strand = location.bs_strand;
}
void match_trace_process(
    matches_t* const matches,
    match_trace_t* const match_trace,
    const locator_t* const locator) {
  // Correct CIGAR (Reverse it if the search was performed in the reverse strand, emulated)
  if (match_trace->emulated_rc_search) {
    matches_cigar_reverse(matches->cigar_vector,
        match_trace->match_alignment.cigar_offset,
        match_trace->match_alignment.cigar_length);
    match_trace->emulated_rc_search = false;
  }
  // Locate-map the match
  match_trace_locate(match_trace,locator);
}
void match_trace_replace(match_trace_t* const match_trace_dst,match_trace_t* const match_trace_src) {
  /* Text (Reference) */
  match_trace_dst->text_trace_offset = match_trace_src->text_trace_offset;
  match_trace_dst->text = match_trace_src->text;
  match_trace_dst->text_length = match_trace_src->text_length;
  /* Match */
  // match_trace_dst->match_trace_offset = match_trace_src->match_trace_offset; // Keep the old
  match_trace_dst->sequence_name = match_trace_src->sequence_name;
  match_trace_dst->strand = match_trace_src->strand;
  match_trace_dst->bs_strand = match_trace_src->bs_strand;
  match_trace_dst->text_position = match_trace_src->text_position;
  match_trace_dst->emulated_rc_search = match_trace_src->emulated_rc_search;
  /* Score */
  match_trace_dst->distance = match_trace_src->distance;
  match_trace_dst->edit_distance = match_trace_src->edit_distance;
  match_trace_dst->swg_score = match_trace_src->swg_score;
  match_trace_dst->mapq_score = match_trace_src->mapq_score;
  /* Alignment */
  match_trace_dst->match_alignment = match_trace_src->match_alignment;
  match_trace_dst->match_scaffold = match_trace_src->match_scaffold;
}
void matches_add_match(
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    const bool preserve_rank,
    select_parameters_t* const select_parameters,
    const alignment_model_t alignment_model,
    match_trace_t** const match_trace_added,
    bool* const match_added,
    bool* const match_replaced,
    mm_stack_t* const mm_stack) {
  // Init
  match_trace->mapq_score = 0;
  // Check duplicates
  uint64_t* dup_match_trace_offset = matches_lookup_match(matches,
      match_trace->match_alignment.match_position,match_trace->match_alignment.effective_length);
  if (dup_match_trace_offset==NULL) {
    // Process match-trace (Correct CIGAR & locate)
    match_trace_process(matches,match_trace,locator);
    // Add the match-trace
    PROF_INC_COUNTER(GP_MATCHES_SE_ADD_NUM_MAPS);
    const uint64_t new_match_trace_offset = vector_get_used(matches->position_matches);
    match_trace_t* new_match_trace;
    vector_alloc_new(matches->position_matches,match_trace_t,new_match_trace);
    *new_match_trace = *match_trace;
    // Index
    matches_index_match(matches,new_match_trace,new_match_trace_offset,mm_stack);
    // Update counters
    matches_counters_add(matches->counters,match_trace->distance,1);
    // Update metrics
    matches_metrics_update(&matches->metrics,match_trace->distance,
        match_trace->edit_distance,match_trace->swg_score); // Update Metrics
    // Preserve rank consistency
    if (preserve_rank) {
      *match_trace_added = matches_sort_preserve_rank_consistency(matches,
          select_parameters,alignment_model,*(new_match_trace->match_trace_offset));
    } else {
      *match_trace_added = new_match_trace;
    }
    // Set added (not replaced)
    *match_added = true;
    *match_replaced = false;
  } else {
    match_trace_t* const dup_match_trace = vector_get_elm(
        matches->position_matches,*dup_match_trace_offset,match_trace_t);
    // Pick up the least distant match
    if (dup_match_trace->distance > match_trace->distance) {
      // Process match-trace (Correct CIGAR & locate)
      match_trace_process(matches,match_trace,locator);
      // Update counters
      matches_counters_sub(matches->counters,dup_match_trace->distance,1);
      matches_counters_add(matches->counters,match_trace->distance,1);
      // Replace the match-trace
      match_trace_replace(dup_match_trace,match_trace);
      // Recompute metrics
      matches_recompute_metrics(matches);
      // Set replaced
      *match_replaced = true;
      // Preserve rank consistency
      if (preserve_rank) {
        *match_trace_added = matches_sort_preserve_rank_consistency(matches,
            select_parameters,alignment_model,*(dup_match_trace->match_trace_offset));
      } else {
        *match_trace_added = dup_match_trace;
      }
    } else {
      // Set not-replaced
      *match_replaced = false;
    }
    // Set not-added (maybe replaced...)
    *match_added = false;
  }
}
bool matches_add_match_trace(
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    mm_stack_t* const mm_stack) {
  match_trace_t* match_trace_added;
  bool match_added, match_replaced;
  matches_add_match(matches,locator,match_trace,false,NULL,
      alignment_model_none,&match_trace_added,&match_added,&match_replaced,mm_stack);
  return match_added;
}
void matches_add_match_trace__preserve_rank(
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    select_parameters_t* const select_parameters,
    const alignment_model_t alignment_model,
    match_trace_t** const match_trace_added,
    bool* const match_added,
    bool* const match_replaced,
    mm_stack_t* const mm_stack) {
  matches_add_match(matches,locator,match_trace,true,select_parameters,
      alignment_model,match_trace_added,match_added,match_replaced,mm_stack);
}
/*
 * Local Matches
 */
void matches_add_local_match_pending(
    matches_t* const matches,
    match_trace_t* const match_trace) {
  match_trace_t* new_match_trace;
  vector_alloc_new(matches->local_matches,match_trace_t,new_match_trace);
  *new_match_trace = *match_trace;
}
void matches_add_pending_local_matches(
    matches_t* const matches,
    const locator_t* const locator,
    mm_stack_t* const mm_stack) {
  const uint64_t num_local_matches = vector_get_used(matches->local_matches);
  match_trace_t* const local_match = vector_get_mem(matches->local_matches,match_trace_t);
  uint64_t i;
  for (i=0;i<num_local_matches;++i) {
    matches_add_match_trace(matches,locator,local_match,mm_stack);
  }
}
/*
 * Matches hints
 */
void matches_hint_allocate_match_trace(
    matches_t* const matches,
    const uint64_t num_matches_trace_to_add) {
  vector_reserve_additional(matches->position_matches,num_matches_trace_to_add);
}
/*
 * Sorting Matches
 */
int match_trace_cmp_distance(const match_trace_t* const a,const match_trace_t* const b) {
  const int distance_diff = (int)a->distance - (int)b->distance;
  if (distance_diff) return distance_diff;
  const int distance_swg = (int)b->swg_score - (int)a->swg_score;
  if (distance_swg) return distance_swg;
  const int distance_edit = (int)a->edit_distance - (int)b->edit_distance;
  if (distance_edit) return distance_edit;
  // Untie using position (helps to stabilize & cmp results)
  return (int)a->match_alignment.match_position - (int)b->match_alignment.match_position;
}
void matches_sort_by_distance(matches_t* const matches) {
  // Sort
  const uint64_t num_matches = vector_get_used(matches->position_matches);
  match_trace_t* const match_trace = vector_get_mem(matches->position_matches,match_trace_t);
  qsort(match_trace,num_matches,sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_distance);
  // Recompute offsets
  uint64_t i;
  for (i=0;i<num_matches;++i) {
    *(match_trace[i].match_trace_offset) = i;
  }
}
int match_trace_cmp_swg_score(const match_trace_t* const a,const match_trace_t* const b) {
  return (int)b->swg_score - (int)a->swg_score;
}
void matches_sort_by_swg_score(matches_t* const matches) {
  // Sort
  const uint64_t num_matches = vector_get_used(matches->position_matches);
  match_trace_t* const match_trace = vector_get_mem(matches->position_matches,match_trace_t);
  qsort(match_trace,num_matches,sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_swg_score);
  // Recompute offsets
  uint64_t i;
  for (i=0;i<num_matches;++i) {
    *(match_trace[i].match_trace_offset) = i;
  }
}
int match_trace_cmp_mapq_score(const match_trace_t* const a,const match_trace_t* const b) {
  return (int)b->mapq_score - (int)a->mapq_score;
}
void matches_sort_by_mapq_score(matches_t* const matches) {
  qsort(vector_get_mem(matches->position_matches,match_trace_t),
      vector_get_used(matches->position_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_mapq_score);
}
int match_trace_cmp_sequence_name__position(const match_trace_t* const a,const match_trace_t* const b) {
  const int cmp = gem_strcmp(a->sequence_name,b->sequence_name);
  return (cmp!=0) ? cmp : ((int)a->text_position - (int)b->text_position);
}
void matches_sort_by_sequence_name__position(matches_t* const matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(matches->position_matches,match_trace_t),
      vector_get_used(matches->position_matches),sizeof(match_trace_t),
      (int (*)(const void *,const void *))match_trace_cmp_sequence_name__position);
}
/*
 * Filters
 */
void matches_filter_by_mapq(matches_t* const matches,const uint8_t mapq_threshold,mm_stack_t* const mm_stack) {
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  match_trace_t* match_in = matches_get_match_trace_buffer(matches);
  match_trace_t* match_out = match_in;
  bool reallocated = false;
  uint64_t i;
  for (i=0;i<num_matches;++i,++match_in) {
    if (match_in->mapq_score >= mapq_threshold) {
      if (match_out != match_in) {
        reallocated = true;
        *match_out = *match_in;
      }
      ++match_out;
    }
  }
  vector_update_used(matches->position_matches,match_out);
  // Because the matches had been reallocated, the indexed-positions are no longer valid
  if (reallocated) {
    matches_index_rebuild(matches,mm_stack);
    matches_recompute_metrics(matches);
  }
}
/*
 * Display
 */
void matches_print(FILE* const stream,matches_t* const matches) {
  tab_fprintf(stream,"[GEM]>Matches\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Class %s\n",matches_class_label[matches->matches_class]);
  tab_fprintf(stream,"=> Counters\t");
  matches_counters_print(stream,matches->counters,matches->max_complete_stratum);
  fprintf(stream,"\n");
  tab_fprintf(stream,"=> Position.Matches %lu\n",vector_get_used(matches->position_matches));
  tab_fprintf(stream,"  => Positions.Hashed.Begin %lu\n",ihash_get_num_elements(matches->begin_pos_matches));
  tab_fprintf(stream,"  => Positions.Hashed.End %lu\n",ihash_get_num_elements(matches->end_pos_matches));
  tab_fprintf(stream,"=> Metrics.SE\n");
  tab_global_inc();
  matches_metrics_print(stream,&matches->metrics);
  tab_global_dec();
  tab_global_dec();
}


