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
  uint64_t num_standard_regions;
  uint64_t total_candidates;
  // Current Region
  uint64_t begin_position;
  // Region Optimization
  uint64_t last_cut;
  uint64_t lo_cut;
  uint64_t hi_cut;
  uint64_t expected_count;
  uint64_t max_steps;
  // Max region length
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
  region_profile->num_standard_regions = 0;
  // Mismatch regions
  region_profile->search_region = mm_stack_calloc(mm_stack,pattern_length,region_search_t,false);
  region_profile->num_search_regions = 0;
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
      (region_profile->filtering_region[0].end == 0 &&
       region_profile->filtering_region[0].start == region_profile->pattern_length &&
       (region_profile->filtering_region[0].hi-region_profile->filtering_region[0].lo) > 0);
}
/*
 * Region Profile Generation
 */
GEM_INLINE void region_profile_generator_close(
    region_profile_generator_t* const rp_generator,
    region_profile_query_t* const rp_query,
    const uint64_t lo,const uint64_t hi) {
  region_search_t* const current_region = rp_generator->region_search + rp_generator->num_regions;
  // Set range
  current_region->end = rp_query->key_position;
  const uint64_t region_length = current_region->start-current_region->end;
  rp_generator->max_region_length = MAX(rp_generator->max_region_length,region_length);
  // Set interval
  current_region->lo = lo;
  current_region->hi = hi;
  const uint64_t candidates_region = hi - lo;
  if (hi - lo <= candidates_region) {
    current_region->type = region_unique;
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
  current_region->start = rp_query->key_position;
  current_region->min = 0;
  rp_generator->last_cut = 0;
  // Region-Query Status
  rp_query->lo = 0;
  rp_query->hi = fm_index_get_length(rp_query->fm_index);
  rank_mquery_new(&rp_query->rank_mquery);
}
GEM_INLINE void region_profile_generator_add_character(
    region_profile_generator_t* const rp_generator,const region_profile_model_t* const profile_model,
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
      region_profile_generator_close(rp_generator,rp_query,lo,hi);
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
      region_profile_generator_close(rp_generator,rp_query,rp_generator->lo_cut,rp_generator->hi_cut);
      region_profile_generator_restart(rp_generator,rp_query);
    }
    return;
  }
  // Zero candidates & (allow zero-regions or no cutting point)
  if (gem_expect_false(allow_zero_regions || rp_generator->last_cut == 0)) {
    region_profile_generator_close(rp_generator,rp_query,lo,hi);
    region_profile_generator_restart(rp_generator,rp_query);
    return;
  }
  // Don't allow zero-regions (restore last cut-point)
  rp_query->key_position = rp_generator->last_cut;
  region_profile_generator_close(rp_generator,rp_query,rp_generator->lo_cut,rp_generator->hi_cut);
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
        (filtering_region[i].start-filtering_region[i].end)<<16 :
        (filtering_region[i].end-filtering_region[i].start);
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
GEM_INLINE void region_profile_fill_gaps(
    region_profile_t* const region_profile,const uint64_t eff_mismatches) {
  GEM_NOT_IMPLEMENTED(); // TODO
//  const bool* const allowed_chars = search_params->allowed_chars;
//  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
//  filtering_region* const filtered_region = region_profile->filtering_region;
//
//  // Prepare memory for the new contiguous regions
//  vector_prepare(mpool->buf3,filtering_region,key_len);
//  filtering_region* c_region = vector_get_mem(mpool->buf3);
//  uint64_t num_regions, last_pos = key_len, next_pos, i;
//  for (i=0,num_regions=0; i<num_filtering_regions; ++num_regions, ++c_region) {
//    next_pos = filtered_region[i].start;
//    if (last_pos == next_pos) { // Contiguous regions
//      c_region->start=filtered_region[i].start; c_region->end=filtered_region[i].end;
//      c_region->type=filtered_region[i].type; c_region->min=filtered_region[i].min;
//      last_pos = filtered_region[i].end;
//      FMI_REGION_SET_NUMBER_BASES(c_region,c_region->start-c_region->end);
//      ++i;
//    } else { // Fill the gap
//      c_region->start=last_pos; c_region->end=next_pos;
//      c_region->type=GAP_REGION; c_region->min=0;
//      last_pos = next_pos;
//      FMI_REGION_PROFILE_COUNT_BASES(key,c_region);
//    }
//  }
//  if (last_pos != 0) {
//    c_region->start=last_pos; c_region->end=0;
//    c_region->type=GAP_REGION; c_region->min=0;
//    FMI_REGION_PROFILE_COUNT_BASES(key,c_region);
//    ++num_regions;
//  }
//
//  // Copy the region profile back
//  c_region = vector_get_mem(mpool->buf3);
//  for (i=0; i<num_regions; ++i) {
//    filtered_region[i].start=c_region[i].start; filtered_region[i].end=c_region[i].end;
//    filtered_region[i].type=c_region[i].type; filtered_region[i].min=c_region[i].min;
//    FMI_REGION_SET_NUMBER_BASES(filtered_region+i,FMI_REGION_GET_NUMBER_BASES(c_region+i))
//  }
//
//  // Substitute filtering region vector
//  region_profile->num_filtering_regions = num_regions;
}
GEM_INLINE void region_profile_extend_last_region(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const bool* const allowed_enc,
    const uint64_t rp_region_type_th) {
  // Tries to extend the last region of the profile to the end of the key (position 0)
  // Merge the tail with the last region
  region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
  if (last_region->end != 0) {
    // Continue the search
    while (last_region->end > 0) {
      if (last_region->hi==last_region->lo) break;
      // Query step
      const uint8_t enc_char = key[--last_region->end];
      if (!allowed_enc[enc_char]) break;
      last_region->lo=bwt_erank(fm_index->bwt,enc_char,last_region->lo);
      last_region->hi=bwt_erank(fm_index->bwt,enc_char,last_region->hi);
    }
    // Extend beyond zero interval
    while (last_region->end > 0 && allowed_enc[key[--last_region->end]]);
    // Adjust the region end (if needed) and type
    const uint64_t count = last_region->hi-last_region->lo;
    if (last_region->type==region_standard && count<=rp_region_type_th) {
      last_region->type = region_unique;
      --(region_profile->num_standard_regions);
    }
  }
}
GEM_INLINE void region_profile_disallow_character(
    region_profile_generator_t* const rp_generator,region_profile_query_t* const rp_query,
    const uint8_t* const key,const bool* const allowed_enc) {
  if (rp_generator->last_cut != 0) {
    ++(rp_query->key_position);
    region_profile_generator_close(rp_generator,rp_query,rp_query->lo,rp_query->hi);
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
  region_profile->max_region_length = rp_generator->max_region_length;
  region_profile->total_candidates = rp_generator->total_candidates;
  if (rp_generator->num_regions == 0) {
    region_search_t* const first_region = rp_generator->region_search;
    if (first_region->start == key_length) {
      first_region->end = 0;
      first_region->lo = rp_query->lo;
      first_region->hi = rp_query->hi;
      region_profile->num_filtering_regions = 1;
      region_profile->num_standard_regions = 1;
    } else {
      region_profile->num_filtering_regions = 0;
      region_profile->num_standard_regions = 0;
    }
  } else {
    // We extend the last region
    region_profile->num_filtering_regions = rp_generator->num_regions;
    region_profile->num_standard_regions = rp_generator->num_standard_regions;
    if (extend_last_region) {
      region_profile_extend_last_region(region_profile,rp_query->fm_index,
          key,allowed_enc,profile_model->region_type_th);
    }
  }
}
/*
 * Region Profile Adaptive
 *
 *   Extracts the adaptive region profile from the given read. That it's, tries
 *   to determine regions of the read which have few matches in the index.
 *
 *   region_th :: Maximum number of matches to consider a region
 *   max_steps :: Maximum number of characters that will be explored trying
 *       to reduce the number of candidates belonging to that region
 *   dec_factor :: Since a region is considered as so, we explore one character a
 *       a time considering if its worth to add that character to that region.
 *       The tradeoff is that we are willing to add that character provided that
 *       the number of candidates gets reduced by 1/rp_dec_factor
 *   region_type_th :: Depending on this threshold there would be 2 kinds of regions
 *       - Regular Regions: Those regions which generates some few results in a exact search (> rp_region_type_th)
 *       - Unique Regions: Those regions which generates very few matches or none (<= rp_region_type_th)
 *   max_regions :: No more than max_regions will be extracted
 *   allow_zero_regions :: Allow a region to have zero candidates
 *
 *   NOTE that if there are no regions this implies that:
 *        (1) There are wildcards which prevents regions generation
 *        (2) There are exact matches
 */
GEM_INLINE void region_profile_generate_adaptive(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,const uint64_t max_regions,
    const bool allow_zero_regions) {
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
      .num_standard_regions = 0,
      .total_candidates = 0,
      .begin_position = 0,
      .last_cut = 0,
      .lo_cut = 0,
      .hi_cut = 0,
      .expected_count = 0,
      .max_steps = 0,
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
      region_profile_disallow_character(&rp_generator,&rp_query,key,allowed_enc);
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
      region_profile_disallow_character(&rp_generator,&rp_query,key,allowed_enc);
      region_length = 0;
      continue;
    }
    // Rank query
    region_profile_generator_query_character(&rp_generator,&rp_query,enc_char);
    ++region_length;
    // Add the character to the region profile
    const uint64_t num_candidates = rp_query.hi-rp_query.lo;
    if (num_candidates <= profile_model->region_th || region_length >= max_region_length) {
      region_profile_generator_close(&rp_generator,&rp_query,rp_query.lo,rp_query.hi);
      region_profile_generator_restart(&rp_generator,&rp_query);
      region_length = 0;
    }
  }
  // Check number of regions
  region_profile_compose(region_profile,&rp_generator,profile_model,
      &rp_query,key,key_length,allowed_enc,false);
  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) { fprintf(stderr,"\n"); }
  PROF_STOP(GP_REGION_PROFILE_ADAPTIVE);
}
GEM_INLINE void region_profile_generate_full_progressive(
    region_profile_t* const region_profile,region_search_t* const base_region,
    const uint64_t start_region,const uint64_t total_regions) {
  region_profile->num_search_regions = total_regions;
  region_profile->search_region = base_region;
  // Sum the degree for the last regions and zero
  uint64_t i, degree_accum;
  for (i=start_region,degree_accum=0; i<total_regions; ++i) {
    degree_accum += base_region[i].degree;
    base_region[i].degree = REGION_FILTER_NONE;
  }
  // Correct the first region degree
  base_region[0].degree += degree_accum;
}
/*
 * Display
 */
GEM_INLINE void region_profile_print(
    FILE* const stream,const region_profile_t* const region_profile,
    const bool sorted,const bool display_misms_regions) {
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
          region->end,region->start,region->hi-region->lo);
    }
  } else {
    REGION_LOCATOR_ITERATE(region_profile,region,position) {
      tab_fprintf(stream,"    [%lu]\ttype=%s\tregion=[%lu,%lu)\tcand=%lu\t\n",
          position,
          region->type==region_unique ? "region_unique" :
              (region->type==region_standard ? "region_standard" : "region_gap"),
          region->end,region->start,region->hi-region->lo);
    }
  }
  if (display_misms_regions) {
    tab_fprintf(stream,"  => Num.Misms.regions %lu\n",region_profile->num_search_regions);
    tab_fprintf(stream,"  => Misms.Regions\n");
    uint64_t i;
    for (i=0;i<region_profile->num_search_regions;++i) {
      const region_search_t* const mismatch_region = region_profile->search_region+i;
      tab_fprintf(stream,"       {%lu}\tloc=[%lu,%lu)\tdegree=%lu\n",
          i,mismatch_region->end,mismatch_region->start,mismatch_region->degree);
    }
  }
  fflush(stream);
}

