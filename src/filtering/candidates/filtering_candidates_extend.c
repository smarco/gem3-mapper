/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
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
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Filtering module provides functions to extend a match as to
 *   find an alignment in the nearby bases
 */

#include "filtering/candidates/filtering_candidates_extend.h"
#include "filtering/candidates/filtering_candidates_align.h"
#include "filtering/region/filtering_region_verify.h"

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
void filtering_candidates_extend_match(
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
  const uint64_t max_extension_error = candidate_pattern->max_extension_error;
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
      candidate_pattern->key_length,max_extension_error,max_template_size,
      &candidate_begin_position,&candidate_end_position);
  const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
  text_trace_t text_trace;
  archive_text_retrieve(
      archive_text,candidate_begin_position,candidate_length,
      false,false,&text_trace,filtering_candidates->mm_allocator);
  PROFILE_STOP(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS,PROFILE_LEVEL);
  /*
   * Verify candidate region (may contain multiple matches)
   */
  const uint64_t candidates_found = filtering_region_verify_extension(
      filtering_candidates,&text_trace,candidate_begin_position,
      candidate_pattern,max_extension_error);
  text_trace_destroy(&text_trace,filtering_candidates->mm_allocator); // Free
  if (candidates_found==0) { PROFILE_STOP(GP_FC_EXTEND_MATCH,PROFILE_LEVEL); return; }
  /*
   * Align accepted candidates
   */
  PROFILE_START(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  matches_t* matches_candidate = (candidate_end==paired_end1) ?
      paired_matches->matches_end1 : paired_matches->matches_end2;
  // Align
  filtering_candidates_align_extended_candidates(
      filtering_candidates,candidate_pattern,matches_candidate);
  PROFILE_STOP(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  PROFILE_STOP(GP_FC_EXTEND_MATCH,PROFILE_LEVEL);
}
