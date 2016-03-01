/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_extend.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates_extend.h"
#include "filtering/filtering_candidates_align.h"
#include "filtering/filtering_region_verify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

void filtering_candidates_compute_extension_region(
    filtering_candidates_t* const filtering_candidates,
    const bool extension_onward,
    const uint64_t extended_begin_position,
    const uint64_t extended_effective_length,
    const uint64_t candidate_key_length,
    const uint64_t max_filtering_error,
    const uint64_t max_template_size,
    uint64_t* const candidate_begin_position,
    uint64_t* const candidate_end_position) {
  locator_t* const locator = filtering_candidates->archive->locator;
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
uint64_t filtering_candidates_extend_match(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const candidate_pattern,
    const match_trace_t* const extended_match,
    paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end,
    mapper_stats_t* const mapper_stats) {
  PROFILE_START(GP_FC_EXTEND_MATCH,PROFILE_LEVEL);
  // Parameters
  archive_t* const archive = filtering_candidates->archive;
  archive_text_t* const archive_text = archive->text;
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  const uint64_t max_filtering_error = candidate_pattern->max_effective_filtering_error;
  /*
   * Retrieve text-candidate
   */
  PROFILE_START(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Compute region on the opposite strand
  const uint64_t extended_effective_length = extended_match->match_alignment.effective_length;
  const uint64_t extended_match_position = archive_text_get_projection(
      archive_text,extended_match->match_alignment.match_position,extended_effective_length);
  const bool search_onward = false;
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  uint64_t candidate_begin_position, candidate_end_position;
  filtering_candidates_compute_extension_region(filtering_candidates,
      search_onward,extended_match_position,extended_effective_length,
      candidate_pattern->key_length,max_filtering_error,max_template_size,
      &candidate_begin_position,&candidate_end_position);
  const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
  const uint64_t text_trace_offset =
      archive_text_retrieve_collection(archive_text,text_collection,
          candidate_begin_position,candidate_length,false,false,mm_stack);
  PROFILE_STOP(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS,PROFILE_LEVEL);
  /*
   * Verify candidate region (may contain multiple matches)
   */
  uint64_t candidates_found = filtering_region_verify_extension(
      filtering_candidates,text_trace_offset,candidate_begin_position,candidate_pattern);
  if (candidates_found==0) { PROFILE_STOP(GP_FC_EXTEND_MATCH,PROFILE_LEVEL); return 0; }
  /*
   * Align accepted candidates
   */
  PROFILE_START(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  matches_t* matches_candidate = (candidate_end==paired_end1) ?
      paired_matches->matches_end1 : paired_matches->matches_end2;
  matches_hint_allocate_match_trace(matches_candidate,candidates_found); // Hint to matches
  // Align
  candidates_found = filtering_candidates_align_candidates(
      filtering_candidates,candidate_pattern,false,true,false,matches_candidate);
  PROFILE_STOP(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  PROFILE_STOP(GP_FC_EXTEND_MATCH,PROFILE_LEVEL);
  // Return number of extended-matches found
  return candidates_found;
}
void filtering_candidates_extend_generate_candidates(
    filtering_candidates_t* const extended_filtering_candidates,
    filtering_candidates_t* const candidate_filtering_candidates,
    const pattern_t* const extended_pattern,
    const pattern_t* const candidate_pattern,
    paired_matches_t* const paired_matches,
    mapper_stats_t* const mapper_stats) {
  // Parameters
  archive_t* const archive = candidate_filtering_candidates->archive;
  archive_text_t* const archive_text = archive->text;
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  // Allocate candidate filtering regions
  const uint64_t num_filtering_regions = vector_get_used(extended_filtering_candidates->filtering_regions);
  filtering_region_t* regions_extended = vector_get_mem(extended_filtering_candidates->filtering_regions,filtering_region_t);
  // Traverse extended-regions and generate candidate-regions
  // Split the extension candidates in chunks and add them to the candidate filtering-regions
  uint64_t n;
  for (n=0;n<num_filtering_regions;++n,++regions_extended) {
    // Compute candidate region boundaries
    const uint64_t extended_eff_begin_position = regions_extended->text_begin_position;
    uint64_t candidate_begin_position, candidate_end_position;
    filtering_candidates_compute_extension_region(
        candidate_filtering_candidates,true,extended_eff_begin_position,
        extended_pattern->key_length+extended_pattern->max_effective_filtering_error,
        candidate_pattern->key_length,candidate_pattern->max_effective_filtering_error,
        max_template_size,&candidate_begin_position,&candidate_end_position);
    const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
    // Add candidate region
    const uint64_t max_error = candidate_pattern->max_effective_filtering_error;
    const uint64_t pattern_length = candidate_pattern->key_length;
    const uint64_t window_length = pattern_length + 2*max_error;
    const uint64_t overlap = pattern_length + max_error;
    const uint64_t begin_position = archive_text_get_projection(archive_text,candidate_begin_position,candidate_length);
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
        regions_candidate.text_begin_position = BOUNDED_SUBTRACTION(current_begin_position,overlap,begin_position);
        regions_candidate.text_end_position = BOUNDED_ADDITION(current_begin_position,overlap,end_position);
        vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
        ++num_chunks_added;
      }
      // Add new chunk
      regions_candidate.text_begin_position = current_begin_position;
      regions_candidate.text_end_position = current_end_position;
      vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
      ++num_chunks_added;
      // Next
      current_begin_position += window_length;
      current_end_position += window_length;
    }
    // Last chunk
    if (num_chunks_added==0 ||
       (current_begin_position < end_position && end_position-current_begin_position > overlap)) {
      regions_candidate.text_begin_position = current_begin_position;
      regions_candidate.text_end_position = end_position;
      vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
    } else if (current_begin_position == begin_position) {
      filtering_region_t* const last_regions_candidate =
          vector_get_last_elm(candidate_filtering_candidates->filtering_regions,filtering_region_t);
      last_regions_candidate->text_end_position = end_position;
    }
  }
}
