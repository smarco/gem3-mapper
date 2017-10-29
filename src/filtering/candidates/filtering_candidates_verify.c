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
 *   Filtering candidates module provides functions to verify filtering-regions
 *   against its corresponding region of text in the index and compute the
 *   distance of the alignment between both
 */

#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/region/filtering_region_verify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Exact-Extend Candidates
 */
void filtering_candidates_extend_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  filtering_region_t** const regions = filtering_candidates_get_regions(filtering_candidates);
  const uint64_t num_regions = filtering_candidates_get_num_regions(filtering_candidates);
  // Extend all discarded regions
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    // Fetch region and retrieve Text
    filtering_region_t* const filtering_region = regions[i];
    filtering_region_retrieve_text(filtering_region,
        pattern,archive_text,filtering_candidates->mm_allocator);
    // Exact extend and region chain scaffold
    text_trace_t* const text_trace = &filtering_region->text_trace;
    if (filtering_region->match_scaffold.scaffold_type <= scaffold_none) {
      match_scaffold_region_chain(
          &filtering_region->match_scaffold,pattern,
          text_trace,matches,filtering_candidates->mm_allocator);
    }
  }
  // Sort by coverage
  filtering_candidates_sort_regions_by_scaffold_coverage(filtering_candidates);
}
void filtering_candidates_extend_discarded_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  filtering_region_t** const regions_discarded = filtering_candidates_get_discarded_regions(filtering_candidates);
  const uint64_t num_regions_discarded = filtering_candidates_get_num_discarded_regions(filtering_candidates);
  // Extend all discarded regions
  uint64_t i;
  for (i=0;i<num_regions_discarded;++i) {
    // Fetch region and retrieve Text
    filtering_region_t* const filtering_region = regions_discarded[i];
    filtering_region_retrieve_text(filtering_region,
        pattern,archive_text,filtering_candidates->mm_allocator);
    // Exact extend and region chain scaffold
    text_trace_t* const text_trace = &filtering_region->text_trace;
    if (filtering_region->match_scaffold.scaffold_type <= scaffold_none) {
      match_scaffold_region_chain(
          &filtering_region->match_scaffold,pattern,
          text_trace,matches,filtering_candidates->mm_allocator);
    }
  }
  // Sort by coverage
  filtering_candidates_sort_discarded_by_scaffold_coverage(filtering_candidates);
}
/*
 * Candidate Verification
 */
uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Check number of filtering regions
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  if (num_filtering_regions==0) return 0;
  PROFILE_START(GP_FC_VERIFICATION,PROFILE_LEVEL);
  PROFILE_START(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const bool candidate_drop_off = search_parameters->candidate_verification.candidate_local_drop_off;
  const verification_strategy_t verification_strategy = search_parameters->candidate_verification.verification_strategy;
  const bool kmer_count_filter = (verification_strategy >= verification_chained);
  // Candidate drop-off filter
  if (candidate_drop_off) {
    filtering_candidates_extend_candidates(filtering_candidates,pattern,matches);
  }
  // Verify all regions
  filtering_region_t** const regions_in = filtering_candidates_get_regions(filtering_candidates);
  filtering_region_t** const regions_out = regions_in;
  filtering_region_t** const regions_discarded =
      filtering_candidates_reserve_discarded_regions(filtering_candidates,num_filtering_regions);
  uint64_t n, num_regions_accepted = 0, num_regions_discarded = 0;
  for (n=0;n<num_filtering_regions;++n) {
    // Verify region
    filtering_region_t* const filtering_region = regions_in[n];
    if (filtering_region_verify(filtering_candidates,filtering_region,kmer_count_filter,pattern)) {
      regions_out[num_regions_accepted++] = filtering_region;
    } else {
      // Heuristic filter
      if (candidate_drop_off && num_regions_accepted==0) {
        const float coverage = filtering_region->match_scaffold.scaffolding_coverage;
        const float key_length = pattern->key_length;
        // Check coverage
        if (coverage/key_length < 1.0) {
          for (;n<num_filtering_regions;++n) { // Discard all regions
            regions_discarded[num_regions_discarded++] = regions_in[n];
          }
          break;
        }
      }
      // Discard
      regions_discarded[num_regions_discarded++] = filtering_region;
    }
  }
  // Update Used
  filtering_candidates_set_num_regions(filtering_candidates,num_regions_accepted);
  filtering_candidates_add_num_discarded_regions(filtering_candidates,num_regions_discarded);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  PROFILE_STOP(GP_FC_VERIFICATION,PROFILE_LEVEL);
  return num_regions_accepted;
}
