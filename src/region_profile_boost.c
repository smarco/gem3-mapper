/*
 * PROJECT: GEMMapper
 * FILE: region_profile_boost.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "region_profile.h"
#include "region_profile_boost.h"
#include "pattern.h"
#include "fm_index_query.h"

/*
 * Debug
 */
#define REGION_PROFILE_DEBUG_PRINT_PROFILE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Region Profile Generator & Query
 */
typedef struct {
  uint64_t key_position;
  uint64_t lo;
  uint64_t hi;
} region_profile_query_t;

/*
 * Region Profile Boost
 */
GEM_INLINE uint64_t region_profile_generate_adaptive_boost_region(
    fm_index_t* const fm_index,const uint8_t* const key,const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,
    const uint64_t region_begin,const uint64_t region_end_max,
    region_search_t* const regions_boosted,uint64_t* const num_regions_boosted,
    uint64_t* const total_candidates,uint64_t* const max_region_length) {
  // Parameters
  const uint64_t region_mayor_th = profile_model->region_th;
  const uint64_t region_minor_th = profile_model->region_type_th;
  const uint64_t fm_index_length = fm_index_get_length(fm_index);
  // Init Query Structures
  region_profile_query_t query;
  query.lo = 0;
  query.hi = fm_index_length;
  query.key_position = region_begin;
  fm_2interval_t fm_2interval;
  fm_2interval.backward_lo = 0;
  fm_2interval.backward_hi = fm_index_length;
  fm_2interval.forward_lo = 0;
  fm_2interval.forward_hi = fm_index_length;
  // Extend-search region
  while (query.key_position < region_end_max) {
    // Get next character
    const uint8_t enc_char = key[query.key_position];
    if (!allowed_enc[enc_char]) break; // Handling wildcards
    query.key_position++;
    // Rank query
    fm_index_2query_forward(fm_index,&fm_2interval,enc_char);
    if (fm_2interval.backward_hi-fm_2interval.backward_lo <= region_minor_th) break;
    query.lo = fm_2interval.backward_lo;
    query.hi = fm_2interval.backward_hi;
  }
  // Add the region to boosted
  const uint64_t num_candidates = fm_2interval.backward_hi-fm_2interval.backward_lo;
  if (num_candidates <= region_mayor_th) {
    region_search_t* const region_boosted = regions_boosted + *num_regions_boosted;
    (*num_regions_boosted)++;
    if (num_candidates == 0 && query.hi-query.lo <= region_mayor_th) {
      region_boosted->begin = region_begin;
      region_boosted->end = query.key_position-1;
      region_boosted->lo = query.lo;
      region_boosted->hi = query.hi;
      region_boosted->type = region_standard;
      *total_candidates += query.hi - query.lo;
    } else {
      region_boosted->begin = region_begin;
      region_boosted->end = query.key_position;
      region_boosted->lo = fm_2interval.backward_lo;
      region_boosted->hi = fm_2interval.backward_hi;
      region_boosted->type = region_standard;
      *total_candidates += num_candidates;
    }
    const uint64_t region_length = region_boosted->end-region_boosted->begin;
    *max_region_length = MAX(*max_region_length,region_length);
    return region_boosted->end;
  } else {
    return region_end_max;
  }
}
GEM_INLINE void region_profile_generate_adaptive_boost(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
  // Parameters
  const double proper_length = fm_index->proper_length;
  const uint64_t min_region_length = (uint64_t)(REGION_MAX_LENGTH_PL_FACTOR*proper_length);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    tab_fprintf(gem_log_get_stream(),"[GEM]>Region.Profile.Generate.Boost\n");
    tab_fprintf(gem_log_get_stream(),"[#%"PRIu64"]",region_profile_num++);
    pattern_enc_print(stderr,key,key_length);
    fprintf(gem_log_get_stream(),"\n");
  }
  // Allocate regions boosted
  mm_stack_push_state(mm_stack);
  region_search_t* const regions_boosted = mm_stack_calloc(mm_stack,key_length,region_search_t,false);
  uint64_t i, num_regions_boosted = 0, total_candidates = 0, max_region_length = 0;
  // Boost cases
  if (region_profile->num_filtering_regions==0) {
    // Boost last region gap
    region_profile_generate_adaptive_boost_region(
        fm_index,key,allowed_enc,profile_model,0,key_length,
        regions_boosted,&num_regions_boosted,&total_candidates,&max_region_length);
  } else {
    // Boost large regions
    for (i=0;i<region_profile->num_filtering_regions;++i) {
      // Check region needs to be boosted
      region_search_t* const region = region_profile->filtering_region + i;
      const uint64_t region_length = region->end-region->begin;
      if (region_length <= min_region_length) {
        max_region_length = MAX(max_region_length,region_length);
        continue; // This region doesn't need to be boosted
      }
      // Boost region
      uint64_t region_end = region_profile_generate_adaptive_boost_region(
          fm_index,key,allowed_enc,profile_model,region->begin,region->end,
          regions_boosted,&num_regions_boosted,&total_candidates,&max_region_length);
      if (region->end-region_end >= proper_length) {
        // (Re)Boost region
        region_profile_generate_adaptive_boost_region(
            fm_index,key,allowed_enc,profile_model,region_end,key_length,
            regions_boosted,&num_regions_boosted,&total_candidates,&max_region_length);
      }
    }
    // Check last-region gap
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    if (last_region->begin >= proper_length) {
      // Boost last region gap
      region_profile_generate_adaptive_boost_region(
          fm_index,key,allowed_enc,profile_model,0,key_length,
          regions_boosted,&num_regions_boosted,&total_candidates,&max_region_length);
    }
  }
  // Compose final region profile
  for (i=0;i<num_regions_boosted;++i) {
    region_profile->filtering_region[i] = regions_boosted[i];
  }
  region_profile->num_filtering_regions = num_regions_boosted;
  region_profile->total_candidates = total_candidates;
  region_profile->num_standard_regions = num_regions_boosted;
  region_profile->num_unique_regions = 0;
  region_profile->num_zero_regions = 0;
  region_profile->max_region_length = max_region_length;
  mm_stack_pop_state(mm_stack,false); // Free
  PROFILE_STOP(GP_REGION_PROFILE_ADAPTIVE,PROFILE_LEVEL);
}
