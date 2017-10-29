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
 *   Region-Profile module provides functions to generate an key partition
 *   based of the entropy of the regions. An adaptive profile determines a
 *   key partition into regions that have few matches in the index.
 *     Note that if the algorithm cannot find any region
 *     could be due to the following reasons
 *       - There are wildcards which prevents region generation
 *       - There are too many exact matches (preventing unique regions)
 */

#include "align/pattern/pattern.h"
#include "filtering/region_profile/region_profile_adaptive.h"
#include "filtering/region_profile/region_profile_schedule.h"

/*
 * Debug
 */
#define REGION_PROFILE_DEBUG_PRINT_PROFILE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Constants
 */
#define REGION_CUTPOINT_NULL UINT64_MAX

/*
 * Region Profile Generation (Query)
 */
void region_profile_generator_save_cut_point(region_profile_generator_t* const generator) {
  generator->last_cut = generator->key_position;
  generator->lo_cut = generator->lo;
  generator->hi_cut = generator->hi;
}
void region_profile_generator_restart(region_profile_generator_t* const generator) {
  region_profile_t* const region_profile = generator->region_profile;
  region_search_t* const current_region = region_profile->filtering_region + region_profile->num_filtering_regions;
  current_region->end = generator->key_position;
  current_region->degree = 0;
  generator->last_cut = REGION_CUTPOINT_NULL;
  // Region-Query Status
  generator->lo = 0;
  generator->hi = fm_index_get_length(generator->fm_index);
  rank_mquery_new(&generator->rank_mquery);
}
/*
 * Region Profile Generation (Region Handling)
 */
void region_profile_generator_init(
    region_profile_generator_t* const generator,
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool allow_zero_regions) {
  // Region profile
  generator->region_profile = region_profile;
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,region_profile->max_expected_regions); // Allocate
  // Region state
  generator->expected_count = 0;
  generator->max_steps = 0;
  // Query
  generator->fm_index = fm_index;
  generator->key = key;
  generator->key_length = key_length;
  generator->allow_zero_regions = allow_zero_regions;
  // Query state
  generator->key_position = key_length;
  region_profile_generator_restart(generator);
}
void region_profile_generator_close_region(
    region_profile_generator_t* const generator,
    const uint64_t key_position,
    const uint64_t lo,
    const uint64_t hi) {
  region_profile_t* const region_profile = generator->region_profile;
  region_search_t* const current_region = region_profile->filtering_region + region_profile->num_filtering_regions;
  // Set range
  current_region->begin = key_position;
  // Set interval/candidates
  current_region->lo = lo;
  current_region->hi = hi;
  ++(region_profile->num_filtering_regions);
}
void region_profile_generator_close_profile(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model) {
  // Check last cut-point
  if (generator->last_cut != REGION_CUTPOINT_NULL) {
    region_profile_generator_close_region(generator,
        generator->last_cut,generator->lo_cut,generator->hi_cut);
  }
  // Check region-profile
  region_profile_t* const region_profile = generator->region_profile;
  if (region_profile->num_filtering_regions == 0) {
    region_search_t* const first_region = region_profile->filtering_region;
    if (first_region->end == generator->key_length) { // Exact Match
      first_region->begin = 0;
      first_region->lo = generator->lo;
      first_region->hi = generator->hi;
      region_profile->num_filtering_regions = 1;
    } else {
      region_profile->num_filtering_regions = 0;
    }
  } else {
    // We extend the last region
    if (generator->allow_zero_regions) {
      region_profile_extend_last_region(
          region_profile,generator->fm_index,generator->key);
    }
  }
  // Select regions & schedule filtering
  region_profile_schedule_exact(region_profile,ALL);
}
bool region_profile_generator_add_character(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model,
    const uint64_t proper_length) {
  // Region lookup status
  const uint64_t lo = generator->lo;
  const uint64_t hi = generator->hi;
  const uint64_t num_candidates = hi-lo;
  // Check number of candidates
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    fprintf(gem_log_get_stream()," %"PRIu64,num_candidates);
  }
  if (num_candidates > profile_model->region_th) return false;
  if (num_candidates > 0) {
    // If we don't have a Cut-Point
    if (generator->last_cut == REGION_CUTPOINT_NULL) {
      region_profile_generator_save_cut_point(generator); // First Cut-Point
      generator->expected_count = num_candidates;
      generator->max_steps = profile_model->max_steps;
      return false;
    }
    // Check Region-Candidates Progress
    generator->expected_count /= profile_model->dec_factor;
    if (num_candidates < generator->expected_count) {
      generator->expected_count = num_candidates; // Dynamic Update
      region_profile_generator_save_cut_point(generator); // Refresh cut point
    }
    // Check maximum steps allowed to optimize region
    --(generator->max_steps);
    if (generator->max_steps == 0) {
      region_profile_generator_close_region(generator,
          generator->last_cut,generator->lo_cut,generator->hi_cut);
      generator->key_position = generator->last_cut; // Restart
      region_profile_generator_restart(generator);
      return true;
    }
    return false;
  } else { // num_candidates == 0
    // Zero candidates & (allow zero-regions or no cutting point)
    if (generator->allow_zero_regions || generator->last_cut == REGION_CUTPOINT_NULL) {
      region_profile_generator_close_region(generator,
          generator->key_position,lo,hi);
    } else {
      // Don't allow zero-regions (or we have a restore last cut-point)
      generator->key_position = generator->last_cut;
      region_profile_generator_close_region(generator,
          generator->last_cut,generator->lo_cut,generator->hi_cut);
    }
    region_profile_generator_restart(generator);
    return true;
  }
}
bool region_profile_generator_disallow_character(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model) {
  bool new_region = false;
  if (generator->last_cut != REGION_CUTPOINT_NULL) {
    region_profile_generator_close_region(generator,
        generator->last_cut,generator->lo_cut,generator->hi_cut);
    new_region = true;
  }
  while (generator->key_position > 0 &&
         generator->key[generator->key_position-1]==ENC_DNA_CHAR_N) {
    --(generator->key_position);
  }
  region_profile_generator_restart(generator);
  return new_region;
}
/*
 * Region Profile Adaptive Iterator
 */
bool region_profile_generator_next_region(
    region_profile_t* const region_profile,
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Parameters
  const uint64_t proper_length = fm_index_get_proper_length(generator->fm_index);
  // Delimit regions
  while (generator->key_position > 0) {
    // Get next character
    --(generator->key_position);
    const uint8_t enc_char = generator->key[generator->key_position];
    // Handling wildcards
    if (enc_char == ENC_DNA_CHAR_N) {
      if (region_profile_generator_disallow_character(generator,profile_model)) {
        PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
        return true; // New region available
      }
    } else {
      // Rank query
      region_profile_query_character(generator->fm_index,
          &generator->rank_mquery,&generator->lo,&generator->hi,enc_char);
      // Add the character to the region profile
      if (region_profile_generator_add_character(generator,profile_model,proper_length)) {
        PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
        return true; // New region available
      }
    }
  }
  // EOI
  region_profile_generator_close_profile(generator,profile_model);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
  return false;
}
/*
 * Region Profile Adaptive
 */
void region_profile_generate_adaptive(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const region_profile_model_t* const profile_model,
    const uint64_t max_regions,
    const bool allow_zero_regions) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    tab_fprintf(gem_log_get_stream(),"[GEM]>Region.Profile.Generate.Adaptive\n");
    tab_fprintf(gem_log_get_stream(),"[#%"PRIu64"]",region_profile_num++);
    pattern_enc_print(stderr,key,key_length);
    fprintf(gem_log_get_stream(),"\n");
    tab_fprintf(gem_log_get_stream(),"[Trace]");
  }
  // Parameters
  const uint64_t proper_length = fm_index_get_proper_length(fm_index);
  const uint64_t max_regions_profiled = MIN(max_regions,region_profile->max_expected_regions);
  // Init
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,
      region_profile,fm_index,key,key_length,allow_zero_regions);
  // Delimit regions
  while (generator.key_position > 0) {
    // Cut-off
    if (generator.region_profile->num_filtering_regions >= max_regions_profiled) {
      PROF_INC_COUNTER(GP_REGION_PROFILE_QUIT_PROFILE);
      break;
    }
    // Get next character
    --(generator.key_position);
    const uint8_t enc_char = key[generator.key_position];
    // Handling wildcards
    if (enc_char == ENC_DNA_CHAR_N) {
      region_profile_generator_disallow_character(&generator,profile_model);
    } else {
      // Rank query
      region_profile_query_character(generator.fm_index,
          &generator.rank_mquery,&generator.lo,&generator.hi,enc_char);
      // Add the character to the region profile
      region_profile_generator_add_character(&generator,profile_model,proper_length);
    }
  }
  region_profile_generator_close_profile(&generator,profile_model);
  // DEBUG
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    fprintf(gem_log_get_stream(),"\n");
  }
}
/*
 * Region Profile Adaptive (limited to extract a minimum number of regions)
 */
void region_profile_generate_adaptive_limited(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const region_profile_model_t* const profile_model) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    tab_fprintf(gem_log_get_stream(),"[GEM]>Region.Profile.Generate.Adaptive\n");
    tab_fprintf(gem_log_get_stream(),"[#%"PRIu64"]",region_profile_num++);
    pattern_enc_print(stderr,key,key_length);
    fprintf(gem_log_get_stream(),"\n");
    tab_fprintf(gem_log_get_stream(),"[Trace]");
  }
  // Parameters
  rank_mtable_t* const rank_mtable = fm_index->rank_table;
  const uint64_t min_matching_depth = rank_mtable->min_matching_depth;
  const uint64_t min_regions = profile_model->num_regions;
  // Init
  const uint64_t max_region_length = key_length/min_regions;
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,region_profile,fm_index,key,key_length,false);
  // Delimit regions
  while (generator.key_position > 0) {
    // Cut-off
    if (generator.region_profile->num_filtering_regions >= region_profile->max_expected_regions) break;
    // Get next character
    --(generator.key_position);
    const uint8_t enc_char = key[generator.key_position];
    // Handling wildcards
    if (enc_char == ENC_DNA_CHAR_N) {
      region_profile_generator_disallow_character(&generator,profile_model);
    } else {
      // Rank query
      region_profile_query_character(generator.fm_index,&generator.rank_mquery,&generator.lo,&generator.hi,enc_char);
      // Add the character to the region profile
      region_search_t* const current_region = region_profile->filtering_region + region_profile->num_filtering_regions;
      const uint64_t num_candidates = generator.hi-generator.lo;
      const uint64_t region_length = current_region->end - generator.key_position;
      if (num_candidates <= profile_model->region_th || region_length >= max_region_length) {
        if (generator.rank_mquery.level < min_matching_depth) {
          rank_mtable_fetch(rank_mtable,&generator.rank_mquery,&generator.lo,&generator.hi);
        }
        region_profile_generator_close_region(&generator,
            generator.key_position,generator.lo,generator.hi);
        region_profile_generator_restart(&generator);
      }
    }
  }
  region_profile_generator_close_profile(&generator,profile_model);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) { fprintf(stderr,"\n"); }
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
