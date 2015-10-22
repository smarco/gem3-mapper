/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_extend.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates_extend.h"
#include "filtering_candidates_align.h"
#include "filtering_region_verify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Constants
 */
#define SPLIT_EXTENSION_WINDOW

GEM_INLINE void filtering_candidates_compute_extension_region(
    const locator_t* const locator,const bool extension_onward,
    const uint64_t extended_begin_position,const uint64_t extended_effective_length,
    const uint64_t candidate_key_length,const uint64_t max_filtering_error,
    const uint64_t max_template_size,uint64_t* const candidate_begin_position,
    uint64_t* const candidate_end_position) {
  locator_interval_t* const locator_interval = locator_lookup_interval(locator,extended_begin_position);
  if (extension_onward) {
    *candidate_begin_position = BOUNDED_SUBTRACTION(extended_begin_position,max_filtering_error,locator_interval->begin_position);
    const uint64_t end_offset = extended_effective_length + max_template_size + candidate_key_length + max_filtering_error;
    *candidate_end_position = BOUNDED_ADDITION(extended_begin_position,end_offset,locator_interval->end_position);
  } else {
    *candidate_end_position = extended_begin_position + extended_effective_length + max_filtering_error;
    if (*candidate_end_position > locator_interval->end_position) {
      *candidate_end_position = locator_interval->end_position;
    }
    const uint64_t begin_offset = max_template_size + candidate_key_length + max_filtering_error;
    *candidate_begin_position = BOUNDED_SUBTRACTION(extended_begin_position,begin_offset,locator_interval->begin_position);
  }
}
/*
 * Pair Extension
 */
GEM_INLINE uint64_t filtering_candidates_extend_match(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,const match_trace_t* const extended_match,
    pattern_t* const candidate_pattern,const as_parameters_t* const candidate_actual_parameters,
    mapper_stats_t* const mapper_stats,paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_FC_EXTEND_MATCH,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = candidate_actual_parameters->search_parameters;
  const uint64_t max_filtering_error = candidate_pattern->max_effective_filtering_error;
  /*
   * Retrieve text-candidate
   */
  PROFILE_START(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Compute region on the opposite strand
  const uint64_t extended_effective_length = extended_match->match_alignment.effective_length;
  const uint64_t extended_match_position = archive_text_get_projection(archive_text,
      extended_match->match_alignment.match_position,extended_effective_length);
  const bool search_onward = false;
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  uint64_t candidate_begin_position, candidate_end_position;
  filtering_candidates_compute_extension_region(
      locator,search_onward,extended_match_position,extended_effective_length,
      candidate_pattern->key_length,max_filtering_error,max_template_size,
      &candidate_begin_position,&candidate_end_position);
  const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
  const uint64_t text_trace_offset = archive_text_retrieve(archive_text,text_collection,
      candidate_begin_position,candidate_length,false,mm_stack);
  PROFILE_STOP(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS,PROFILE_LEVEL);
  /*
   * Verify candidate region (may contain multiple matches)
   */
  PROFILE_START(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  uint64_t candidates_found = filtering_region_verify_extension(
      filtering_candidates->filtering_regions,filtering_candidates->verified_regions,
      text_collection,text_trace_offset,candidate_begin_position,search_parameters,candidate_pattern);
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATE_LENGTH,candidate_length);
  PROFILE_STOP(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  if (candidates_found==0) { PROFILE_STOP(GP_FC_EXTEND_MATCH,PROFILE_LEVEL); return 0; }
  /*
   * Align accepted candidates
   */
  PROFILE_START(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  matches_t* matches_candidate = (candidate_end==paired_end1) ? paired_matches->matches_end1 : paired_matches->matches_end2;
  matches_hint_allocate_match_trace(matches_candidate,candidates_found); // Hint to matches
  // Align
  candidates_found = filtering_candidates_align_accepted_regions(filtering_candidates,
      archive_text,locator,text_collection,candidate_pattern,false,
      candidate_actual_parameters,false,matches_candidate,mm_stack);
  PROFILE_STOP(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  PROFILE_STOP(GP_FC_EXTEND_MATCH,PROFILE_LEVEL);
  // Return number of extended-matches found
  return candidates_found;
}
#ifdef SPLIT_EXTENSION_WINDOW
GEM_INLINE void filtering_candidates_extend_generate_candidates(
    filtering_candidates_t* const extended_filtering_candidates,
    filtering_candidates_t* const candidate_filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    const pattern_t* const extended_pattern,const pattern_t* const candidate_pattern,
    const search_parameters_t* const search_parameters,mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches,mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  // Allocate candidate filtering regions
  const uint64_t num_filtering_regions = vector_get_used(extended_filtering_candidates->filtering_regions);
  filtering_region_t* regions_extended = vector_get_mem(extended_filtering_candidates->filtering_regions,filtering_region_t);
  // Traverse extended-regions and generate candidate-regions
  // Split the extension candidates in chunks and add them to the candidate filtering-regions
  uint64_t n;
  for (n=0;n<num_filtering_regions;++n,++regions_extended) {
    // Compute candidate region boundaries
    const uint64_t extended_eff_begin_position = regions_extended->begin_position;
    uint64_t candidate_begin_position, candidate_end_position;
    filtering_candidates_compute_extension_region(
        archive->locator,true,extended_eff_begin_position,
        extended_pattern->key_length+extended_pattern->max_effective_filtering_error,
        candidate_pattern->key_length,candidate_pattern->max_effective_filtering_error,
        max_template_size,&candidate_begin_position,&candidate_end_position);
    const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
    // Add candidate region
    const uint64_t max_error = candidate_pattern->max_effective_filtering_error;
    const uint64_t pattern_length = candidate_pattern->key_length;
    const uint64_t window_length = pattern_length + 2*max_error;
    const uint64_t overlap = pattern_length + max_error;
    const uint64_t begin_position = archive_text_get_projection(archive->text,candidate_begin_position,candidate_length);
    const uint64_t end_position = begin_position + candidate_length;
    uint64_t current_begin_position = begin_position, current_end_position = begin_position + window_length;
    // Init template
    filtering_region_t regions_candidate;
    regions_candidate.status = filtering_region_unverified; // Newly created region (unverified)
    match_scaffold_init(&regions_candidate.match_scaffold);
    uint64_t num_chunks_added = 0;
    while (current_end_position <= end_position) {
      // Add overlapping candidate
      if (current_begin_position != begin_position) {
        regions_candidate.begin_position = BOUNDED_SUBTRACTION(current_begin_position,overlap,begin_position);
        regions_candidate.end_position = BOUNDED_ADDITION(current_begin_position,overlap,end_position);
        vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
        ++num_chunks_added;
      }
      // Add new chunk
      regions_candidate.begin_position = current_begin_position;
      regions_candidate.end_position = current_end_position;
      vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
      ++num_chunks_added;
      // Next
      current_begin_position += window_length;
      current_end_position += window_length;
    }
    // Last chunk
    if (num_chunks_added==0 || (current_begin_position < end_position && end_position-current_begin_position > overlap)) {
      regions_candidate.begin_position = current_begin_position;
      regions_candidate.end_position = end_position;
      vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
    } else if (current_begin_position == begin_position) {
      filtering_region_t* const last_regions_candidate =
          vector_get_last_elm(candidate_filtering_candidates->filtering_regions,filtering_region_t);
      last_regions_candidate->end_position = end_position;
    }
  }
}
#else
GEM_INLINE void filtering_candidates_process_extension_candidates(
    filtering_candidates_t* const extended_filtering_candidates,
    filtering_candidates_t* const candidate_filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,const pattern_t* const extended_pattern,
    const pattern_t* const candidate_pattern,const search_parameters_t* const search_parameters,
    paired_matches_t* const paired_matches,mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t min_template_size = mapper_stats_template_length_get_expected_min(mapper_stats);
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  // Allocate candidate filtering regions
  const uint64_t num_filtering_regions = vector_get_used(extended_filtering_candidates->filtering_regions);
  filtering_region_t* regions_extended = vector_get_mem(extended_filtering_candidates->filtering_regions,filtering_region_t);
  vector_reserve_additional(candidate_filtering_candidates->filtering_regions,num_filtering_regions);
  filtering_region_t* regions_candidate = vector_get_mem(candidate_filtering_candidates->filtering_regions,filtering_region_t);
  // Traverse extended-regions and generate candidate-regions
  uint64_t n;
  for (n=0;n<num_filtering_regions;++n,++regions_extended,++regions_candidate) {
    // Compute candidate region boundaries
    const uint64_t extended_eff_begin_position = regions_extended->begin_position;
    uint64_t candidate_begin_position, candidate_end_position;
    filtering_candidates_compute_extension_region(
        archive->locator,true,extended_eff_begin_position,
        extended_pattern->key_length+extended_pattern->max_effective_filtering_error,
        candidate_pattern->key_length,candidate_pattern->max_effective_filtering_error,
        min_template_size,max_template_size,&candidate_begin_position,&candidate_end_position);
    const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
    // Add candidate region
    regions_candidate->status = filtering_region_unverified; // Newly created region (unverified)
    match_scaffold_init(&regions_candidate->match_scaffold);
    regions_candidate->begin_position =
        archive_text_get_projection(archive->text,candidate_begin_position,candidate_length);
    regions_candidate->end_position = regions_candidate->begin_position+candidate_length;
    regions_candidate->num_regions_matching = 0;
    regions_candidate->regions_matching = NULL;
    regions_candidate->coverage = 0;
  }
  vector_add_used(candidate_filtering_candidates->filtering_regions,num_filtering_regions);
}
#endif
