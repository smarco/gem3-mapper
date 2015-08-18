/*
 * PROJECT: GEMMapper
 * FILE: region_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "region_profile.h"
#include "mapper_profile.h"
#include "pattern.h"

/*
 * Debug
 */
#define REGION_PROFILE_DEBUG_PRINT_PROFILE GEM_DEEP_DEBUG

/*
 * Region Profile Generator & Query
 */
typedef struct {
  // Region Profile
  region_search_t* region_search;
  uint64_t num_regions;
  // Current Region
  uint64_t begin_position;
  // Region Optimization
  uint64_t last_cut;
  uint64_t lo_cut;
  uint64_t hi_cut;
  uint64_t expected_count;
  uint64_t max_steps;
  // Stats
  uint64_t num_standard_regions;
  uint64_t num_unique_regions;
  uint64_t num_zero_regions;
  uint64_t total_candidates;
  uint64_t max_region_length;
} region_profile_generator_t;
typedef struct {
  fm_index_t* const fm_index;
  uint64_t key_position;
  rank_mquery_t rank_mquery;
  uint64_t lo;
  uint64_t hi;
} region_profile_query_t;

/*
 * Setup
 */
GEM_INLINE void region_profile_new(
    region_profile_t* const region_profile,const uint64_t pattern_length,
    mm_stack_t* const mm_stack) {
  // Allocate memory for the region profile
  region_profile->pattern_length = pattern_length;
  // Filtering regions
  region_profile->filtering_region = mm_stack_calloc(mm_stack,pattern_length,region_search_t,false);
  region_profile->num_filtering_regions = 0;
  // Region Partition Properties
  region_profile->errors_allowed = 0;
  // Locator for region sorting
  region_profile->loc = mm_stack_calloc(mm_stack,pattern_length,region_locator_t,false);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t region_get_num_regions(region_profile_t* const region_profile) {
  // All regions (region_unique + region_standard + region_gap)
  return region_profile->num_filtering_regions;
}
GEM_INLINE bool region_profile_has_exact_matches(region_profile_t* const region_profile) {
  // Condition (sufficient but not necessary) to determine
  //   if there are exact matches (taken into account the first region)
  return (region_profile->num_filtering_regions==1) &&
      (region_profile->filtering_region[0].begin == 0 &&
       region_profile->filtering_region[0].end == region_profile->pattern_length &&
       (region_profile->filtering_region[0].hi-region_profile->filtering_region[0].lo) > 0);
}
/*
 * Region Profile Generation
 */
GEM_INLINE void region_profile_generator_close(
    region_profile_generator_t* const rp_generator,region_profile_query_t* const rp_query,
    const region_profile_model_t* const profile_model,const uint64_t lo,const uint64_t hi) {
  region_search_t* const current_region = rp_generator->region_search + rp_generator->num_regions;
  // Set range
  current_region->begin = rp_query->key_position;
  const uint64_t region_length = current_region->end-current_region->begin;
  rp_generator->max_region_length = MAX(rp_generator->max_region_length,region_length);
  // Set interval
  current_region->lo = lo;
  current_region->hi = hi;
  const uint64_t candidates_region = hi - lo;
  if (candidates_region <= profile_model->region_type_th) {
    current_region->type = region_unique;
    if (candidates_region==0) ++(rp_generator->num_zero_regions);
  } else {
    current_region->type = region_standard;
    ++(rp_generator->num_standard_regions);
  }
  // Set candidates
  rp_generator->total_candidates += candidates_region;
  ++(rp_generator->num_regions);
}
GEM_INLINE void region_profile_generator_save_cut_point(
    region_profile_generator_t* const rp_generator,region_profile_query_t* const rp_query) {
  rp_generator->last_cut = rp_query->key_position;
  rp_generator->lo_cut = rp_query->lo;
  rp_generator->hi_cut = rp_query->hi;
}
GEM_INLINE void region_profile_generator_restart(
    region_profile_generator_t* const rp_generator,region_profile_query_t* const rp_query) {
  region_search_t* const current_region = rp_generator->region_search + rp_generator->num_regions;
  current_region->end = rp_query->key_position;
  current_region->min = 0;
  rp_generator->last_cut = 0;
  // Region-Query Status
  rp_query->lo = 0;
  rp_query->hi = fm_index_get_length(rp_query->fm_index);
  rank_mquery_new(&rp_query->rank_mquery);
}
GEM_INLINE void region_profile_generator_add_character(
    region_profile_generator_t* const rp_generator,
    const region_profile_model_t* const profile_model,
    region_profile_query_t* const rp_query,const bool allow_zero_regions) {
  // Region lookup status
  const uint64_t lo = rp_query->lo;
  const uint64_t hi = rp_query->hi;
  const uint64_t key_position = rp_query->key_position;
  const uint64_t num_candidates = hi-lo;
  // Check number of candidates
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) { fprintf(stderr," %lu",num_candidates); }
  if (gem_expect_true(num_candidates > profile_model->region_th)) return;
  if (num_candidates > 0) {
    // End of the read reached
    if (gem_expect_false(key_position == 0)) {
      region_profile_generator_close(rp_generator,rp_query,profile_model,lo,hi);
      region_profile_generator_restart(rp_generator,rp_query);
      return;
    }
    // If we don't have a Cut-Point
    if (rp_generator->last_cut == 0) {
      region_profile_generator_save_cut_point(rp_generator,rp_query); // First Cut-Point
      rp_generator->expected_count = num_candidates;
      rp_generator->max_steps = profile_model->max_steps;
      return;
    }
    // Check Region-Candidates Progress
    rp_generator->expected_count /= profile_model->dec_factor;
    if (num_candidates <= rp_generator->expected_count || num_candidates <= profile_model->region_type_th) {
      region_profile_generator_save_cut_point(rp_generator,rp_query); // Refresh cut point
    }
    // Check maximum steps allowed to optimize region
    --(rp_generator->max_steps);
    if (rp_generator->max_steps == 0) {
      rp_query->key_position = rp_generator->last_cut;
      region_profile_generator_close(rp_generator,rp_query,profile_model,rp_generator->lo_cut,rp_generator->hi_cut);
      region_profile_generator_restart(rp_generator,rp_query);
    }
    return;
  }
  // Zero candidates & (allow zero-regions or no cutting point)
  if (gem_expect_false(allow_zero_regions || rp_generator->last_cut == 0)) {
    region_profile_generator_close(rp_generator,rp_query,profile_model,lo,hi);
    region_profile_generator_restart(rp_generator,rp_query);
    return;
  }
  // Don't allow zero-regions (restore last cut-point)
  rp_query->key_position = rp_generator->last_cut;
  region_profile_generator_close(rp_generator,rp_query,profile_model,rp_generator->lo_cut,rp_generator->hi_cut);
  region_profile_generator_restart(rp_generator,rp_query);
}
GEM_INLINE void region_profile_generator_query_character(
    region_profile_generator_t* const rp_generator,
    region_profile_query_t* const rp_query,const uint8_t enc_char) {
  if (!rank_mquery_is_exhausted(&rp_query->rank_mquery)) {
    rank_mquery_t* const rank_mquery = &rp_query->rank_mquery;
    rank_mtable_t* const rank_mtable = rp_query->fm_index->rank_table;
    const uint64_t min_matching_depth = rank_mtable->min_matching_depth;
    rank_mquery_add_char(rank_mtable,rank_mquery,enc_char);
    if (rank_mquery->level >= min_matching_depth) {
      rank_mtable_fetch(rank_mtable,rank_mquery,&rp_query->lo,&rp_query->hi);
    }
  } else {
    bwt_t* const bwt = rp_query->fm_index->bwt;
    const uint64_t lo = rp_query->lo;
    const uint64_t hi = rp_query->hi;
    if (gem_expect_false(bwt_is_same_bucket(lo,hi))) {
      bwt_erank_interval(bwt,enc_char,lo,hi,&rp_query->lo,&rp_query->hi);
    } else {
      rp_query->lo = bwt_erank(bwt,enc_char,lo);
      rp_query->hi = bwt_erank(bwt,enc_char,hi);
    }
  }
}
/*
 * Utils
 */
int region_profile_locator_cmp(const region_locator_t* const a,const region_locator_t* const b) {
  return (int)a->value - (int)b->value;
}
GEM_INLINE void region_profile_locator_sort(region_locator_t* const loc,const uint64_t num_regions) {
  qsort(loc,num_regions,sizeof(region_locator_t),(int (*)(const void *,const void *))region_profile_locator_cmp);
}
GEM_INLINE void region_profile_sort_by_estimated_mappability(region_profile_t* const region_profile) {
  // Sort the regions w.r.t to the number of candidates
  const uint64_t num_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  region_locator_t* const loc = region_profile->loc;
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    loc[i].id = i;
    loc[i].value = (filtering_region[i].type == region_standard) ?
        (filtering_region[i].end-filtering_region[i].begin)<<16 :
        (filtering_region[i].hi-filtering_region[i].lo);
  }
  region_profile_locator_sort(loc,num_regions);
}
GEM_INLINE void region_profile_sort_by_candidates(region_profile_t* const region_profile) {
  // Sort the regions w.r.t to the number of candidates
  const uint64_t num_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  region_locator_t* const loc = region_profile->loc;
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    loc[i].id = i;
    loc[i].value = filtering_region[i].hi-filtering_region[i].lo;
  }
  region_profile_locator_sort(loc,num_regions);
}
GEM_INLINE void region_profile_extend_last_region(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const bool* const allowed_enc,
    const uint64_t rp_region_type_th) {
  // Tries to extend the last region of the profile to the end of the key (position 0)
  // Merge the tail with the last region
  region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
  if (last_region->begin != 0) {
    // Continue the search
    while (last_region->begin > 0) {
      if (last_region->hi==last_region->lo) break;
      // Query step
      const uint8_t enc_char = key[--last_region->begin];
      if (!allowed_enc[enc_char]) break;
      last_region->lo=bwt_erank(fm_index->bwt,enc_char,last_region->lo);
      last_region->hi=bwt_erank(fm_index->bwt,enc_char,last_region->hi);
    }
    // Extend beyond zero interval
    while (last_region->begin > 0 && allowed_enc[key[--last_region->begin]]);
    // Adjust the region end (if needed) and type
    const uint64_t count = last_region->hi-last_region->lo;
    if (last_region->type==region_standard && count<=rp_region_type_th) {
      last_region->type = region_unique;
      --(region_profile->num_standard_regions);
    }
  }
}
GEM_INLINE void region_profile_disallow_character(
    region_profile_generator_t* const rp_generator,
    region_profile_query_t* const rp_query,
    const region_profile_model_t* const profile_model,
    const uint8_t* const key,const bool* const allowed_enc) {
  if (rp_generator->last_cut != 0) {
    ++(rp_query->key_position);
    region_profile_generator_close(rp_generator,rp_query,profile_model,rp_query->lo,rp_query->hi);
    --(rp_query->key_position);
  }
  while (rp_query->key_position > 0 && !allowed_enc[key[rp_query->key_position-1]]) {
    --(rp_query->key_position);
  }
  region_profile_generator_restart(rp_generator,rp_query);
}
GEM_INLINE void region_profile_compose(
    region_profile_t* const region_profile,region_profile_generator_t* const rp_generator,
    const region_profile_model_t* const profile_model,region_profile_query_t* const rp_query,
    const uint8_t* const key,const uint64_t key_length,
    const bool* const allowed_enc,const bool extend_last_region) {
  if (rp_generator->num_regions == 0) {
    region_search_t* const first_region = rp_generator->region_search;
    if (first_region->end == key_length) { // Exact Match
      first_region->begin = 0;
      first_region->lo = rp_query->lo;
      first_region->hi = rp_query->hi;
      region_profile->num_filtering_regions = 1;
      region_profile->num_standard_regions = 1;
      region_profile->num_unique_regions = 0;
      region_profile->num_zero_regions = 0;
    } else {
      region_profile->num_filtering_regions = 0;
      region_profile->num_standard_regions = 0;
      region_profile->num_unique_regions = 0;
      region_profile->num_zero_regions = 0;
    }
  } else {
    region_profile->num_filtering_regions = rp_generator->num_regions;
    region_profile->num_standard_regions = rp_generator->num_standard_regions;
    region_profile->num_unique_regions = rp_generator->num_unique_regions;
    region_profile->num_zero_regions = rp_generator->num_zero_regions;
    // Add information about the last region
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    rp_generator->max_region_length = MAX(rp_generator->max_region_length,last_region->begin);
    if (extend_last_region) {
      // We extend the last region
      region_profile_extend_last_region(region_profile,rp_query->fm_index,
          key,allowed_enc,profile_model->region_type_th);
    }
  }
  region_profile->max_region_length = rp_generator->max_region_length;
  region_profile->total_candidates = rp_generator->total_candidates;
}
/*
 * Region Profile Adaptive
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
GEM_INLINE void region_profile_generate_adaptive(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const uint64_t key_length,
    const bool* const allowed_enc,const region_profile_model_t* const profile_model,
    const uint64_t max_regions,const bool allow_zero_regions) {
  PROF_START(GP_REGION_PROFILE_ADAPTIVE);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    fprintf(stderr,"[%lu]",region_profile_num++);
    uint64_t i;
    for (i=0;i<key_length;++i) fprintf(stderr,"%c",dna_decode(key[i]));
    fprintf(stderr,"\n");
  }
  // Init
  region_profile_generator_t rp_generator = {
      .region_search = region_profile->filtering_region,
      .num_regions = 0,
      .total_candidates = 0,
      .begin_position = 0,
      .last_cut = 0,
      .lo_cut = 0,
      .hi_cut = 0,
      .expected_count = 0,
      .max_steps = 0,
      .num_standard_regions = 0,
      .num_unique_regions = 0,
      .num_zero_regions = 0,
      .total_candidates = 0,
      .max_region_length = 0,
  };
  region_profile_query_t rp_query = {
      .fm_index = fm_index,
      .key_position = key_length,
      .lo = 0,
      .hi = 0,
  };
  rank_mquery_new(&rp_query.rank_mquery);
  // Delimit regions
  region_profile_generator_restart(&rp_generator,&rp_query);
  while (rp_query.key_position > 0) {
    // Cut-off
    if (rp_generator.num_regions >= max_regions) { PROF_INC_COUNTER(GP_REGION_PROFILE_QUIT_PROFILE); break; }
    // Get next character
    --(rp_query.key_position);
    const uint8_t enc_char = key[rp_query.key_position];
    // Handling wildcards
    if (!allowed_enc[enc_char]) {
      region_profile_disallow_character(&rp_generator,&rp_query,profile_model,key,allowed_enc);
      continue;
    }
    // Rank query
    region_profile_generator_query_character(&rp_generator,&rp_query,enc_char);
    // Add the character to the region profile
    region_profile_generator_add_character(&rp_generator,profile_model,&rp_query,allow_zero_regions);
  }
  // Check number of regions
  region_profile_compose(region_profile,&rp_generator,profile_model,
      &rp_query,key,key_length,allowed_enc,allow_zero_regions);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) { fprintf(stderr,"\n"); }
  PROF_STOP(GP_REGION_PROFILE_ADAPTIVE);
}
GEM_INLINE void region_profile_generate_adaptive_limited(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,const uint64_t min_regions) {
  PROF_START(GP_REGION_PROFILE_ADAPTIVE);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    fprintf(stderr,"[%lu]",region_profile_num++);
    pattern_enc_print(stderr,key,key_length);
    fprintf(stderr,"\n");
  }
  // Init
  const uint64_t max_region_length = key_length/min_regions;
  region_profile_generator_t rp_generator = {region_profile->filtering_region, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  region_profile_query_t rp_query = {
      .fm_index = fm_index,
      .key_position = key_length,
      .lo = 0,
      .hi = 0,
  };
  rank_mquery_new(&rp_query.rank_mquery);
  // Delimit regions
  uint64_t region_length = 0;
  region_profile_generator_restart(&rp_generator,&rp_query);
  while (rp_query.key_position > 0) {
    // Get next character
    --(rp_query.key_position);
    const uint8_t enc_char = key[rp_query.key_position];
    // Handling wildcards
    if (!allowed_enc[enc_char]) {
      region_profile_disallow_character(&rp_generator,&rp_query,profile_model,key,allowed_enc);
      region_length = 0;
      continue;
    }
    // Rank query
    region_profile_generator_query_character(&rp_generator,&rp_query,enc_char);
    ++region_length;
    // Add the character to the region profile
    const uint64_t num_candidates = rp_query.hi-rp_query.lo;
    if (num_candidates <= profile_model->region_th || region_length >= max_region_length) {
      region_profile_generator_close(&rp_generator,&rp_query,profile_model,rp_query.lo,rp_query.hi);
      region_profile_generator_restart(&rp_generator,&rp_query);
      region_length = 0;
    }
  }
  region_profile_compose(region_profile,&rp_generator,profile_model,&rp_query,key,key_length,allowed_enc,false);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) { fprintf(stderr,"\n"); }
  PROF_STOP(GP_REGION_PROFILE_ADAPTIVE);
}
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
  region_profile_query_t rp_query;
  rp_query.lo = 0;
  rp_query.hi = fm_index_length;
  rp_query.key_position = region_begin;
  fm_2interval_t fm_2interval;
  fm_2interval.backward_lo = 0;
  fm_2interval.backward_hi = fm_index_length;
  fm_2interval.forward_lo = 0;
  fm_2interval.forward_hi = fm_index_length;
  // Extend-search region
  while (rp_query.key_position < region_end_max) {
    // Get next character
    const uint8_t enc_char = key[rp_query.key_position];
    if (!allowed_enc[enc_char]) break; // Handling wildcards
    rp_query.key_position++;
    // Rank query
    fm_index_2query_forward(fm_index,&fm_2interval,enc_char);
    if (fm_2interval.backward_hi-fm_2interval.backward_lo <= region_minor_th) break;
    rp_query.lo = fm_2interval.backward_lo;
    rp_query.hi = fm_2interval.backward_hi;
  }
  // Add the region to boosted
  const uint64_t num_candidates = fm_2interval.backward_hi-fm_2interval.backward_lo;
  if (num_candidates <= region_mayor_th) {
    region_search_t* const region_boosted = regions_boosted + *num_regions_boosted;
    (*num_regions_boosted)++;
    if (num_candidates == 0 && rp_query.hi-rp_query.lo <= region_mayor_th) {
      region_boosted->begin = region_begin;
      region_boosted->end = rp_query.key_position-1;
      region_boosted->lo = rp_query.lo;
      region_boosted->hi = rp_query.hi;
      region_boosted->type = region_standard;
      *total_candidates += rp_query.hi - rp_query.lo;
    } else {
      region_boosted->begin = region_begin;
      region_boosted->end = rp_query.key_position;
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
  PROF_START(GP_REGION_PROFILE_ADAPTIVE);
  // Parameters
  const double proper_length = fm_index->proper_length;
  const uint64_t min_region_length = (uint64_t)(REGION_MAX_LENGTH_PL_FACTOR*proper_length);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    fprintf(stderr,"[%lu]",region_profile_num++);
    pattern_enc_print(stderr,key,key_length);
    fprintf(stderr,"\n");
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
  PROF_STOP(GP_REGION_PROFILE_ADAPTIVE);
}
/*
 * Display
 */
GEM_INLINE void region_profile_print(
    FILE* const stream,const region_profile_t* const region_profile,const bool sorted) {
  tab_fprintf(stream,"[GEM]>Region.Profile\n");
  tab_fprintf(stream,"  => Pattern.length %lu\n",region_profile->pattern_length);
  tab_fprintf(stream,"  => Num.Filtering.Regions %lu\n",region_profile->num_filtering_regions);
  tab_fprintf(stream,"  => Num.Standard.Regions %lu\n",region_profile->num_standard_regions);
  tab_fprintf(stream,"  => Errors.allowed %lu\n",region_profile->errors_allowed);
  tab_fprintf(stream,"  => Filtering.Regions\n");
  if (!sorted) {
    REGION_PROFILE_ITERATE(region_profile,region,position) {
      tab_fprintf(stream,"    [%lu]\ttype=%s\tregion=[%lu,%lu)\tcand=%lu\n",
          position,
          region->type==region_unique ? "region_unique" :
              (region->type==region_standard ? "region_standard" : "region_gap"),
          region->begin,region->end,region->hi-region->lo);
    }
  } else {
    REGION_LOCATOR_ITERATE(region_profile,region,position) {
      tab_fprintf(stream,"    [%lu]\ttype=%s\tregion=[%lu,%lu)\tcand=%lu\t\n",
          position,
          region->type==region_unique ? "region_unique" :
              (region->type==region_standard ? "region_standard" : "region_gap"),
          region->begin,region->end,region->hi-region->lo);
    }
  }
  fflush(stream);
}

