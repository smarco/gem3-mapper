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
 *   Approximate-String-Matching (ASM) module encapsulating
 *   the basic stages in a step-wise approach. Mainly used to
 *   horizontally process a batch of searches (doing each stage
 *   for all searches at the same time before progressing to
 *   the next search stage)
 */

#include "approximate_search/approximate_search_stages.h"
#include "approximate_search/approximate_search_stepwise.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_verify_candidates.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "filtering/candidates/filtering_candidates_align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Region Profile Utils
 */
void approximate_search_stepwise_region_profile_adaptive_compute(
    approximate_search_t* const search) {
  PROF_START(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL);
  // Re-Compute region profile
  search->processing_state = asearch_processing_state_begin;
  approximate_search_region_profile(search);
  if (search->processing_state==asearch_processing_state_no_regions) return;
  // Set State
  search->processing_state = asearch_processing_state_region_profiled;
  PROF_STOP(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL);
}
/*
 * AM Stepwise :: Region Profile
 */
void approximate_search_stepwise_region_profile_static_generate(
    approximate_search_t* const search) {
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search Start. Check basic cases
        approximate_search_begin(search);
        break;
      case asearch_stage_filtering_adaptive:
        approximate_search_region_profile_static_partition(search);
        return;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
void approximate_search_stepwise_region_profile_static_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  if (search->processing_state == asearch_processing_state_region_partitioned) {
    approximate_search_region_profile_static_buffered_copy(search,gpu_buffer_fmi_ssearch);
  }
}
void approximate_search_stepwise_region_profile_static_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  // Retrieve regions-interval (compute if not enabled)
  if (search->processing_state == asearch_processing_state_region_partitioned) {
    if (gpu_buffer_fmi_ssearch->fmi_search_enabled) {
      approximate_search_region_profile_static_buffered_retrieve(search,gpu_buffer_fmi_ssearch);
    } else {
      approximate_search_region_profile_static_compute(search);
    }
  }
  // Compute Region Profile Adaptively (unsuccessful cases)
  if (search->processing_state == asearch_processing_state_no_regions) {
    approximate_search_stepwise_region_profile_adaptive_compute(search);
    if (search->processing_state == asearch_processing_state_no_regions) return; // Corner cases
  }
  // Check exact matches & limit the number of matches
  if (region_profile_has_exact_matches(&search->region_profile)) {
    approximate_search_generate_candidates_limit_exact_matches(search);
  }
}
void approximate_search_stepwise_region_profile_adaptive_generate(
    approximate_search_t* const search) {
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search Start. Check basic cases
        approximate_search_begin(search);
        break;
      case asearch_stage_filtering_adaptive:
        search->search_stage = asearch_stage_filtering_adaptive;
        return;
      case asearch_stage_end:
    	return;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
void approximate_search_stepwise_region_profile_adaptive_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  if (search->search_stage==asearch_stage_filtering_adaptive) {
    approximate_search_region_profile_adaptive_buffered_copy(search,gpu_buffer_fmi_asearch);
  }
}
void approximate_search_stepwise_region_profile_adaptive_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  // Retrieve regions-interval (compute if not enabled)
  if (search->search_stage==asearch_stage_filtering_adaptive) {
    if (gpu_buffer_fmi_asearch->fmi_search_enabled) {
      approximate_search_region_profile_adaptive_buffered_retrieve(search,gpu_buffer_fmi_asearch);
    } else {
      approximate_search_stepwise_region_profile_adaptive_compute(search);
    }
    // Check exact matches & limit the number of matches
    if (region_profile_has_exact_matches(&search->region_profile)) {
      approximate_search_generate_candidates_limit_exact_matches(search);
    }
    // Set state profiled
    search->processing_state = asearch_processing_state_region_profiled;
  }
}
/*
 * AM Stepwise :: Decode Candidates
 */
void approximate_search_stepwise_decode_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  if (search->processing_state == asearch_processing_state_region_profiled) {
    approximate_search_generate_candidates_buffered_copy(search,gpu_buffer_fmi_decode);
  }
}
void approximate_search_stepwise_decode_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  if (search->processing_state == asearch_processing_state_region_profiled) {
    approximate_search_generate_candidates_buffered_retrieve(search,gpu_buffer_fmi_decode);
  }
}
/*
 * AM Stepwise :: Kmer-filter
 */
void approximate_search_stepwise_kmer_filter_copy(
    approximate_search_t* const search,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_kmer_filter_buffered_copy(search,gpu_buffer_kmer_filter);
  }
}
void approximate_search_stepwise_kmer_filter_retrieve(
    approximate_search_t* const search,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  if (search->search_stage==asearch_stage_filtering_adaptive) {
    if (search->processing_state == asearch_processing_state_candidates_processed) {
      approximate_search_kmer_filter_buffered_retrieve(search,gpu_buffer_kmer_filter);
    }
  }
}
/*
 * AM Stepwise :: BPM-Distance
 */
void approximate_search_stepwise_bpm_distance_copy(
    approximate_search_t* const search,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_bpm_distance_buffered_copy(search,gpu_buffer_bpm_distance);
  }
}
void approximate_search_stepwise_bpm_distance_retrieve(
    approximate_search_t* const search,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  if (search->search_stage==asearch_stage_filtering_adaptive) {
    if (search->processing_state == asearch_processing_state_candidates_processed) {
      approximate_search_bpm_distance_buffered_retrieve(search,gpu_buffer_bpm_distance);
    }
  }
}
/*
 * AM Stepwise :: BPM-Align
 */
void approximate_search_stepwise_bpm_align_update(
    approximate_search_t* const search) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_bpm_align_buffered_update(search);
  }
}
void approximate_search_stepwise_bpm_align_copy(
    approximate_search_t* const search,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_bpm_align_buffered_copy(search,gpu_buffer_bpm_align);
  }
}
void approximate_search_stepwise_bpm_align_retrieve(
    approximate_search_t* const search,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    matches_t* const matches) {
  if (search->search_stage==asearch_stage_filtering_adaptive) {
    if (search->processing_state == asearch_processing_state_candidates_processed) {
      approximate_search_bpm_align_buffered_retrieve(search,gpu_buffer_bpm_align,matches);
    }
    search->search_stage = asearch_stage_filtering_adaptive_finished;
  }
}
