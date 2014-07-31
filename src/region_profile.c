/*
 * PROJECT: GEMMapper
 * FILE: region_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "region_profile.h"

/*
 * Setup
 */
GEM_INLINE void region_profile_new(
    region_profile_t* const region_profile,pattern_t* const pattern,
    mm_stack_t* const mm_stack) {
  // Allocate memory for the region profile
  const uint64_t max_length = pattern->key_length;
  // Pattern
  region_profile->pattern_length = max_length;
  // Filtering regions
  region_profile->filtering_region = mm_stack_calloc(mm_stack,max_length,filtering_region_t,false);
  region_profile->num_filtering_regions = 0;
  region_profile->num_standard_regions = 0;
  // Mismatch regions
  region_profile->mismatch_region = mm_stack_calloc(mm_stack,max_length,mismatch_region_t,false);
  region_profile->num_mismatch_region = 0;
  // Region Partition Properties
  region_profile->misms_required = 0;
  // Locator for region sorting
  region_profile->loc = mm_stack_calloc(mm_stack,max_length,region_locator_t,false);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t region_get_num_regions(region_profile_t* const region_profile) {
  // TODO
  return 0;
}
GEM_INLINE bool region_profile_has_exact_matches(region_profile_t* const region_profile) {
  // Condition to determine if there are exact matches (taken into account the first region)
  return (region_profile->num_filtering_regions>0) &&
      (region_profile->filtering_region[0].end == 0 &&
       region_profile->filtering_region[0].start == region_profile->pattern_length &&
       (region_profile->filtering_region[0].hi-region_profile->filtering_region[0].lo) > 0);
}
/*
 * Region Profile Generation (Fixed)
 *
 *   Extracts a region profile ...
 */
GEM_INLINE void region_profile_generate_fixed(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,pattern_t* const pattern,
    bool* const allowed_enc,const uint64_t num_regions) {
  // TODO
}
/*
 * Region Profile Generation (Adaptive)
 *
 *   Extracts the adaptive region profile from the given read.
 *   That it's, tries to determine regions of the read which have few matches in the index.
 *
 *   rp_region_th :: Maximum number of matches to consider a region
 *   rp_max_steps :: Maximum number of characters that will be explored trying
 *       to reduce the number of candidates belonging to that region
 *   rp_dec_factor :: Since a region is considered as so, we explore one character a
 *       a time considering if its worth to add that character to that region.
 *       The tradeoff is that we are willing to add that character provided that
 *       the number of candidates gets reduced by 1/rp_dec_factor
 *   rp_region_type_th :: Depending on this threshold there would be 2 kinds of regions
 *       - Zero Regions: Those regions which generates some few results in a exact search (> rp_region_type_th)
 *       - Non-Zero Regions: Those regions which generates very few matches or none (<= rp_region_type_th)
 *   max_regions :: No more than max_regions will be extracted
 *
 *   NOTE that if there are no regions this implies that:
 *        (1) There are wildcards which prevents regions generation
 *        (2) There are exact matches
 *      Also note that:
 *        (region_profile->num_filtering_regions==0 && wildcards==0) => ExactMatches
 *      But it could be that:
 *        (region_profile->num_filtering_regions>0) AND ExactMatches
 */
#define REGION_CLOSE(pos,loPos,hiPos) { \
  regions[num_regions].end=pos; \
  regions[num_regions].lo=loPos; \
  regions[num_regions].hi=hiPos; \
  if (hiPos-loPos<=rp_region_type_th) { \
    regions[num_regions].type = region_unique; \
  } else { \
    regions[num_regions].type = region_standard; \
    ++num_standard_regions; \
  } \
  ++num_regions; \
  if (num_regions >= max_regions) { \
    PROF_INC_COUNTER(GP_QUIT_PROFILE); \
    break;\
  } \
}
#define REGION_RESTART(restart_position) { \
  regions[num_regions].start=restart_position; \
  regions[num_regions].min=REGION_FILTER_NONE; \
  lo=0; hi=bwt_length; \
  rank_mquery_new(&rank_mquery); \
  last_cut=0; \
}
#define REGION_SAVE_CUT_POINT() last_cut = key_len; hi_cut = hi; lo_cut = lo
GEM_INLINE void region_profile_generate_adaptive(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,pattern_t* const pattern,const bool* const allowed_enc,
    const uint64_t rp_region_th,const uint64_t rp_max_steps,
    const uint64_t rp_dec_factor,const uint64_t rp_region_type_th,const uint64_t max_regions) {
  PROF_BEGIN(GP_REGION_PROFILE);
  filtering_region_t* const regions = region_profile->filtering_region;
  uint64_t num_regions = 0, num_standard_regions = 0, last_cut = 0;
  uint64_t lo, hi, hi_cut, lo_cut, expected_count, max_steps;

  rank_mquery_t rank_mquery;
  rank_mquery_new(&rank_mquery);

  // Initial region
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  uint8_t* const key = pattern->key;
  uint64_t key_len = pattern->key_length;
  REGION_RESTART(key_len);
  while (key_len > 0) {
    --key_len;

    // Handling wildcards
    if (!allowed_enc[key[key_len]]) {
      if (last_cut != 0) REGION_CLOSE(key_len+1,lo,hi);
      while (key_len > 0 && !allowed_enc[key[key_len-1]]) --key_len;
      REGION_RESTART(key_len);
      continue;
    }

    // Rank query
    const uint8_t enc_char = key[key_len];
    if (rank_mquery_is_exhausted(&rank_mquery)
        /* && lookupLevel > MIN_LOOKUP_FETCH(properLength/4) */) {
      rank_mquery_add_char(&rank_mquery,enc_char);
    } else {
      lo = bwt_erank(fm_index->bwt,enc_char,lo);
      hi = bwt_erank(fm_index->bwt,enc_char,hi);
    }

    // Determine region scope
    const uint64_t num_candidates = hi-lo;
    if (gem_expect_true(num_candidates > rp_region_th)) continue;

    // We have candidates
    if (num_candidates > 0) {
      if (last_cut == 0) {
        REGION_SAVE_CUT_POINT();
        expected_count = num_candidates; max_steps = rp_max_steps;
      } else {
        expected_count /= rp_dec_factor;
        if (num_candidates <= expected_count || num_candidates<=rp_region_type_th) { // Store cut point
          REGION_SAVE_CUT_POINT();
        }
        if (--max_steps == 0) {
          key_len = last_cut;
          REGION_CLOSE(key_len,lo_cut,hi_cut);
          REGION_RESTART(key_len);
        } else if (key_len == 0) {
          REGION_CLOSE(key_len, lo, hi);
          REGION_RESTART(key_len);
        }
      }
    } else {
      REGION_CLOSE(key_len, lo, hi);
      REGION_RESTART(key_len);
    }
  }

  // Store number of regions
  region_profile->num_filtering_regions = num_regions;
  region_profile->num_standard_regions = num_standard_regions;

  // We leave a ghost region
  regions[num_regions].end=0;
  regions[num_regions].lo=lo;
  regions[num_regions].hi=hi;

  PROF_END(GP_REGION_PROFILE);
}
GEM_INLINE void region_profile_generate_full_progressive(
    region_profile_t* const region_profile,mismatch_region_t* const base_region,
    const uint64_t start_region,const uint64_t total_regions) {
  region_profile->num_mismatch_region = total_regions;
  region_profile->mismatch_region = base_region;
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
 * Region Profile Scheduling
 */
GEM_INLINE void region_profile_schedule_exact_filtering(
    region_profile_t* const region_profile,const uint64_t regions_required) {
  /*
   * PRE: (region_profile->num_filtering_regions >= regions_required)
   * Schedules the @regions_required regions with fewer number of candidates
   *   to be filtered up to zero mismatches.
   */
  const uint64_t num_regions = region_profile->num_filtering_regions;
  filtering_region_t* const filtering_region = region_profile->filtering_region;
  region_locator_t* const loc = region_profile->loc;
  // Check the number of regions in the profile
  int64_t i;
  if (num_regions > regions_required) {
    /*
     * More regions than required (2 strategies)
     *   (*) - Filter only those required which have the less number of candidates
     *       - Filter all of them, but only check those
     *           candidates which more than 1 region (H-samples)  // TODO
     */
    // Sort by number of candidates
    region_profile_sort_by_candidates(region_profile);
    // Filter up to 0-mismatches those regions with the less number of candidates
    for (i=num_regions-1;i>=(num_regions-regions_required);--i) {
      filtering_region[loc[i].id].min = REGION_FILTER_DEGREE_ZERO;
    }
    for (;i>=0;--i) {
      filtering_region[loc[i].id].min = REGION_FILTER_NONE;
    }
  } else {
    for (i=0;i<num_regions;++i) {
      loc[i].id = i; // Fill locator to guide the filtering
      filtering_region[i].min = REGION_FILTER_DEGREE_ZERO;
    }
  }
}
/*
 * Region Profile Utils
 */
GEM_INLINE void region_profile_locator_sort(region_locator_t* const loc,const uint64_t num_regions) {
  // Very simple sort function for region locator based on its value
  int64_t i, j;
  region_locator_t elm;
  for (i=1; i<num_regions; ++i) {
    elm = loc[i]; j=i-1;
    while (j>=0 && (elm.value > loc[j].value)) {
      loc[j+1] = loc[j]; --j;
    }
    if (j<i-1) loc[j+1] = elm;
  }
}
GEM_INLINE void region_profile_sort_by_estimated_mappability(region_profile_t* const region_profile) {
  // Sort the regions w.r.t to the number of candidates
  const uint64_t num_regions = region_profile->num_filtering_regions;
  filtering_region_t* const filtering_region = region_profile->filtering_region;
  region_locator_t* const loc = region_profile->loc;
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    loc[i].id = i;
    loc[i].value = (filtering_region[i].type == region_standard) ?
        filtering_region[i].start-filtering_region[i].end : 0;
  }
  region_profile_locator_sort(loc,num_regions);
}
GEM_INLINE void region_profile_sort_by_candidates(region_profile_t* const region_profile) {
  // Sort the regions w.r.t to the number of candidates
  const uint64_t num_regions = region_profile->num_filtering_regions;
  filtering_region_t* const filtering_region = region_profile->filtering_region;
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
/*
 *
 */
GEM_INLINE void region_profile_extend_first_region(
    region_profile_t* const region_profile,fm_index_t* const fm_index,const uint64_t rp_region_type_th) {

}
/*
 * PRE: region_profile->num_regions > 0
 * Tries to extend the last region of the profile to the end of the key (position 0)
 */
GEM_INLINE void region_profile_extend_last_region(
    region_profile_t* const region_profile,fm_index_t* const fm_index,const uint64_t rp_region_type_th) {

//// Merge the tail with the last region
//filtering_region* const last_region = region_profile->filtering_region +
//    (region_profile->num_filtering_regions-1);
//if (last_region->end != 0) {
//  // Continue the search
//  last_region->end = fmi_bsearch_continue(
//      a,key,key_len,allowed_repl,last_region->hi,last_region->lo,
//      last_region->end,0,&(last_region->hi),&(last_region->lo));
//  // Adjust the region end (if needed) and type
//  const uint64_t count = last_region->hi-last_region->lo;
//  if (count==0) last_region->end = 0;
//  if (last_region->type==ZERO_REGION && count<=rp_region_type_th) {
//    last_region->type = NON_ZERO_REGION;
//    --(region_profile->num_zero_filtering_regions);
//  }
//}
}
/*
 * Display
 */
GEM_INLINE void region_profile_print(FILE* const stream,const region_profile_t* const region_profile) {

//#ifdef GEM_SHOW_PROFILES
//#define FMI_SHOW_REGION_PROFILE(REGION_PROFILE, REMARKS...) {
//  filtering_region* ppp_region = REGION_PROFILE->filtering_region;
//  int64_t iii, mmm = REGION_PROFILE->num_filtering_regions;
//  fprintf(stderr, REMARKS);
//  if (mmm==0)  fprintf(stderr, " {no regions}");
//  for (iii=mmm-1; iii>=0; --iii) {
//    fprintf(stderr, "[%lu,%lu]{min=%lu/max=%lu} (degree=%lu;cand=%lu)   <|>   ",
//        ppp_region[iii].end, ppp_region[iii].start,
//        ppp_region[iii].min, ppp_region[iii].max, ppp_region[iii].degree,
//        ppp_region[iii].hi-ppp_region[iii].lo);
//  }
//  fprintf(stderr, "\n");
//}
//#define FMI_SHOW_MISMS_PROFILE(REGION_PROFILE, REMARKS...) {
//  mismatch_region* ppp_region = REGION_PROFILE->mismatch_region;
//  int64_t iii, mmm = REGION_PROFILE->num_mismatch_region;
//  fprintf(stderr, REMARKS);
//  if (mmm==0)  fprintf(stderr, " {no regions}");
//  for (iii=mmm-1; iii>=0; --iii) {
//    fprintf(stderr, "[%lu,%lu](d=%lu)   <-->   ",
//        ppp_region[iii].end, ppp_region[iii].start, ppp_region[iii].degree);
//  }
//  fprintf(stderr, "\n");
//}
}
/*
 * Stats
 */
GEM_INLINE region_profile_stats_t* region_profile_stats_new() {
  // TODO
  return NULL;
}
GEM_INLINE void region_profile_stats_delete(region_profile_stats_t* const region_profile_stats) {
  // TODO
}
GEM_INLINE void region_profile_stats_record(region_profile_stats_t* const region_profile_stats,region_profile_t* const region_profile) {
  // TODO
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////




/*
 * Assigns to @filtered_region profile a continuous set of regions,
 * handling wildcards and regions in-between filtered regions
 */
#define FMI_REGION_GET_NUMBER_BASES(region) (region)->hi
#define FMI_REGION_SET_NUMBER_BASES(region,number_bases) (region)->hi = number_bases;
#define FMI_REGION_PROFILE_COUNT_BASES(key,region) { \
  uint64_t z, count = 0; \
  for (z=region->end; z<region->start; z++)  { \
    if (allowed_chars[key[z]]) ++count; \
  } \
  FMI_REGION_SET_NUMBER_BASES(region,count); \
}

#define FMI_ADD_MISMS_REGION(i,last_pos,size_chunk) \
  misms_region[i].start = last_pos; \
  last_pos -= size_chunk; \
  misms_region[i].end = last_pos






