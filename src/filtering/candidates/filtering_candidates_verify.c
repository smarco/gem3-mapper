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
 * Candidate Verification
 */
uint64_t filtering_candidates_verify_filtering_regions(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern) {
  PROFILE_START(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  // Traverse all regions (text-space)
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  filtering_region_t** const regions_in = filtering_candidates_get_regions(filtering_candidates);
  filtering_region_t** const regions_out = regions_in;
  filtering_region_t** const regions_discarded =
      filtering_candidates_reserve_discarded_regions(filtering_candidates,num_filtering_regions);
  uint64_t n, num_regions_accepted = 0, num_regions_out = 0, num_regions_discarded = 0;
  for (n=0;n<num_filtering_regions;++n) {
    // Check region status (Skip other than unverified)
    filtering_region_t* const filtering_region = regions_in[n];
    if (filtering_region->status!=filtering_region_unverified) {
      regions_out[num_regions_out++] = filtering_region;
    } else {
      // Verify region
      if (filtering_region_verify(filtering_candidates,filtering_region,pattern,true)) {
        regions_out[num_regions_out++] = filtering_region;
        ++num_regions_accepted;
      } else {
        regions_discarded[num_regions_discarded++] = filtering_region;
      }
    }
  }
  // Update Used
  filtering_candidates_set_num_regions(filtering_candidates,num_regions_out);
  filtering_candidates_add_num_discarded_regions(filtering_candidates,num_regions_discarded);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  return num_regions_accepted;
}
uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern) {
  // Check number of filtering regions
  uint64_t pending_candidates = filtering_candidates_get_num_regions(filtering_candidates);
  if (pending_candidates==0) return 0;
  PROFILE_START(GP_FC_VERIFICATION,PROFILE_LEVEL);
  // Verify candidates
  pending_candidates = filtering_candidates_verify_filtering_regions(filtering_candidates,pattern);
  PROFILE_STOP(GP_FC_VERIFICATION,PROFILE_LEVEL);
  return pending_candidates;
}
