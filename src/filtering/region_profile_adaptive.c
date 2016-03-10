/*
 * PROJECT: GEMMapper
 * FILE: region_profile_adaptive.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *
 *   Region Profile Adaptive
 *
 *   Extracts the adaptive region profile from the given read.
 *   Roughly speaking, tries to determine regions of the read which have
 *   few matches in the index. Note that if the algorithm cannot find any region
 *   could be due to the following reasons
 *     - There are wildcards which prevents regions generation
 *     - There are too many exact matches (preventing unique regions)
 *
 *   region_th  - Maximum number of matches allow to determine a region
 *   max_steps  - Maximum number of characters that will be explored after
 *                reaching @region_th trying to reduce the number of candidates
 *                of the region
 *   dec_factor - Once the number of candidates of the region is below @region_th,
 *                the algorithm will expand the region by one character as long as the
 *                total number of candidates of that region is reduced by a factor of @dec_factor
 *   region_type_th - Depending on the number of candidates of the region we classify them into
 *                    regular regions and unique regions.
 *     Regular Regions - Regions with more candidates than @rp_region_type_th (> rp_region_type_th)
 *     Unique Regions  - Regions with less candidates than @rp_region_type_th (<= rp_region_type_th)
 *   max_regions - No more than max_regions will be generated
 *   allow_zero_regions - Allow a region to have zero candidates
 *
 */

#include "filtering/region_profile_adaptive.h"
#include "data_structures/pattern.h"

/*
 * Debug
 */
#define REGION_PROFILE_DEBUG_PRINT_PROFILE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

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
  generator->region_length = 0;
  generator->last_cut = 0;
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
    const bool* const allowed_enc,
    const bool allow_zero_regions) {
  // Region Profile
  generator->region_profile = region_profile;
  region_profile_clear(region_profile);
  // Region State
  generator->region_length = 0;
  generator->last_cut = 0;
  generator->lo_cut = 0;
  generator->hi_cut = 0;
  generator->expected_count = 0;
  generator->max_steps = 0;
  // Query
  generator->fm_index = fm_index;
  generator->key = key;
  generator->key_length = key_length;
  generator->allowed_enc = allowed_enc;
  generator->allow_zero_regions = allow_zero_regions;
  // Query state
  generator->key_position = key_length;
  region_profile_generator_restart(generator);
  // Mappabiliy
  generator->mappability_p_acc = 0.0;
  generator->mappability_p_samples = 0;
  generator->mappability_2p_acc = 0.0;
  generator->mappability_2p_samples = 0;
}
void region_profile_generator_close_region(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model,
    const uint64_t lo,
    const uint64_t hi) {
  region_profile_t* const region_profile = generator->region_profile;
  region_search_t* const current_region = region_profile->filtering_region + region_profile->num_filtering_regions;
  // Set range
  current_region->begin = generator->key_position;
  region_profile->max_region_length = MAX(region_profile->max_region_length,generator->region_length);
  // Set interval
  current_region->lo = lo;
  current_region->hi = hi;
  const uint64_t candidates_region = hi - lo;
  if (candidates_region <= profile_model->region_type_th) {
    current_region->type = region_unique;
    if (candidates_region==0) ++(region_profile->num_zero_regions);
  } else {
    current_region->type = region_standard;
    ++(region_profile->num_standard_regions);
  }
  // Set candidates
  region_profile->total_candidates += candidates_region;
  ++(region_profile->num_filtering_regions);
}
void region_profile_generator_close_profile(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model) {
  region_profile_t* const region_profile = generator->region_profile;
  if (region_profile->num_filtering_regions == 0) {
    region_search_t* const first_region = region_profile->filtering_region;
    if (first_region->end == generator->key_length) { // Exact Match
      first_region->begin = 0;
      first_region->lo = generator->lo;
      first_region->hi = generator->hi;
      region_profile->num_filtering_regions = 1;
      region_profile->num_standard_regions = 1;
      region_profile->num_unique_regions = 0;
      region_profile->num_zero_regions = 0;
      region_profile->total_candidates = generator->hi - generator->lo;
    } else {
      region_profile->num_filtering_regions = 0;
      region_profile->num_standard_regions = 0;
      region_profile->num_unique_regions = 0;
      region_profile->num_zero_regions = 0;
      region_profile->total_candidates = 0;
    }
  } else {
    // We extend the last region
    if (generator->allow_zero_regions) {
      region_profile_extend_last_region(region_profile,generator->fm_index,
          generator->key,generator->allowed_enc,profile_model->region_type_th);
    }
    // Add information about the last region
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    region_profile->max_region_length = MAX(region_profile->max_region_length,last_region->begin);
  }
  // Compute Mappability
  if (generator->mappability_p_samples > 0) {
    region_profile->mappability_p = generator->mappability_p_acc / (double)(2*generator->mappability_p_samples);
  }
  if (generator->mappability_2p_samples > 0) {
    region_profile->mappability_2p = generator->mappability_2p_acc / (double)(2*generator->mappability_2p_samples);
  }
}
bool region_profile_generator_add_character(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model,
    const uint64_t proper_length) {
  // Region lookup status
  const uint64_t lo = generator->lo;
  const uint64_t hi = generator->hi;
  const uint64_t num_candidates = hi-lo;
  // Record Mappability
  ++(generator->region_length);
  if (generator->region_length == proper_length) {
    if (num_candidates > 0) generator->mappability_p_acc += log2((double)num_candidates); // (x/2.0)
    ++(generator->mappability_p_samples);
  } else if (generator->region_length == 2*proper_length) {
    if (num_candidates > 0) generator->mappability_2p_acc += log2((double)num_candidates); // (x/2.0)
    ++(generator->mappability_2p_samples);
  }
  // Check number of candidates
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    fprintf(gem_log_get_stream()," %"PRIu64,num_candidates);
  }
  if (num_candidates > profile_model->region_th) return false;
  if (num_candidates > 0) {
    // End of the read reached
    if (gem_expect_false(generator->key_position == 0)) {
      region_profile_generator_close_region(generator,profile_model,lo,hi);
      region_profile_generator_restart(generator);
      return true;
    }
    // If we don't have a Cut-Point
    if (generator->last_cut == 0) {
      region_profile_generator_save_cut_point(generator); // First Cut-Point
      generator->expected_count = num_candidates;
      generator->max_steps = profile_model->max_steps;
      return false;
    }
    // Check Region-Candidates Progress
    generator->expected_count /= profile_model->dec_factor;
    if (num_candidates <= generator->expected_count || num_candidates <= profile_model->region_type_th) {
      region_profile_generator_save_cut_point(generator); // Refresh cut point
    }
    // Check maximum steps allowed to optimize region
    --(generator->max_steps);
    if (generator->max_steps == 0) {
      generator->key_position = generator->last_cut;
      region_profile_generator_close_region(generator,profile_model,generator->lo_cut,generator->hi_cut);
      region_profile_generator_restart(generator);
      return true;
    }
    return false;
  } else { // num_candidates == 0
    // Zero candidates & (allow zero-regions or no cutting point)
    if (gem_expect_false(generator->allow_zero_regions || generator->last_cut == 0)) {
      region_profile_generator_close_region(generator,profile_model,lo,hi);
      region_profile_generator_restart(generator);
      return true;
    }
    // Don't allow zero-regions (restore last cut-point)
    generator->key_position = generator->last_cut;
    region_profile_generator_close_region(generator,profile_model,generator->lo_cut,generator->hi_cut);
    region_profile_generator_restart(generator);
    return true;
  }
}
bool region_profile_generator_disallow_character(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model) {
  bool new_region = false;
  if (generator->last_cut != 0) {
    ++(generator->key_position);
    region_profile_generator_close_region(generator,profile_model,generator->lo,generator->hi);
    --(generator->key_position);
    new_region = true;
  }
  while (generator->key_position > 0 && !generator->allowed_enc[generator->key[generator->key_position-1]]) {
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
  PROFILE_START(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
  // Parameters
  const uint64_t proper_length = fm_index_get_proper_length(generator->fm_index);
  // Delimit regions
  while (generator->key_position > 0) {
    // Get next character
    --(generator->key_position);
    const uint8_t enc_char = generator->key[generator->key_position];
    // Handling wildcards
    if (!generator->allowed_enc[enc_char]) {
      if (region_profile_generator_disallow_character(generator,profile_model)) {
        PROFILE_STOP(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
        return true; // New region available
      }
    } else {
      // Rank query
      region_profile_query_character(generator->fm_index,
          &generator->rank_mquery,&generator->lo,&generator->hi,enc_char);
      // Add the character to the region profile
      if (region_profile_generator_add_character(generator,profile_model,proper_length)) {
        PROFILE_STOP(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
        return true; // New region available
      }
    }
  }
  // EOI
  region_profile_generator_close_profile(generator,profile_model);
  PROFILE_STOP(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
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
    const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,
    const uint64_t max_regions,
    const bool allow_zero_regions) {
  PROFILE_START(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
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
  // Init
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,region_profile,fm_index,key,key_length,allowed_enc,allow_zero_regions);
  // Delimit regions
  while (generator.key_position > 0) {
    // Cut-off
    if (generator.region_profile->num_filtering_regions >= max_regions) {
      PROF_INC_COUNTER(GP_REGION_PROFILE_QUIT_PROFILE);
      break;
    }
    // Get next character
    --(generator.key_position);
    const uint8_t enc_char = key[generator.key_position];
    // Handling wildcards
    if (!allowed_enc[enc_char]) {
      region_profile_generator_disallow_character(&generator,profile_model);
    } else {
      // Rank query
      region_profile_query_character(generator.fm_index,&generator.rank_mquery,&generator.lo,&generator.hi,enc_char);
      // Add the character to the region profile
      region_profile_generator_add_character(&generator,profile_model,proper_length);
    }
  }
  region_profile_generator_close_profile(&generator,profile_model);
  // DEBUG
  PROFILE_STOP(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
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
    const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,
    const uint64_t min_regions) {
  PROFILE_START(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    tab_fprintf(gem_log_get_stream(),"[GEM]>Region.Profile.Generate.Adaptive\n");
    tab_fprintf(gem_log_get_stream(),"[#%"PRIu64"]",region_profile_num++);
    pattern_enc_print(stderr,key,key_length);
    fprintf(gem_log_get_stream(),"\n");
    tab_fprintf(gem_log_get_stream(),"[Trace]");
  }
  // Init
  const uint64_t max_region_length = key_length/min_regions;
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,region_profile,fm_index,key,key_length,allowed_enc,true);
  // Delimit regions
  region_profile_generator_restart(&generator);
  while (generator.key_position > 0) {
    // Get next character
    --(generator.key_position);
    const uint8_t enc_char = key[generator.key_position];
    // Handling wildcards
    if (!allowed_enc[enc_char]) {
      region_profile_generator_disallow_character(&generator,profile_model);
    } else {
      // Rank query
      region_profile_query_character(generator.fm_index,&generator.rank_mquery,&generator.lo,&generator.hi,enc_char);
      // Add the character to the region profile
      const uint64_t num_candidates = generator.hi-generator.lo;
      if (num_candidates <= profile_model->region_th || generator.region_length >= max_region_length) {
        region_profile_generator_close_region(&generator,profile_model,generator.lo,generator.hi);
        region_profile_generator_restart(&generator);
      }
    }
  }
  region_profile_generator_close_profile(&generator,profile_model);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) { fprintf(stderr,"\n"); }
  PROFILE_STOP(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
}
