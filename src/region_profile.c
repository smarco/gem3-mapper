/*
 * PROJECT: GEMMapper
 * FILE: region_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "region_profile.h"
#include "mapper_profile.h"

/*
 * Debug
 */
#define REGION_PROFILE_DEBUG_PRINT_PROFILE false

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
 * Tries to extend the last region of the profile to the end of the key (position 0)
 */
GEM_INLINE void region_profile_extend_last_region(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const bool* const allowed_enc,
    const uint64_t rp_region_type_th) {
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
/*
 * Region Profile Generation
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
    PROF_INC_COUNTER(GP_REGION_PROFILE_QUIT_PROFILE); \
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
#define REGION_SAVE_CUT_POINT(pos,hi,lo) last_cut = pos; hi_cut = hi; lo_cut = lo
///*
// *
// */
//GEM_INLINE void region_profile_generate_overlapped(
//    region_profile_t* const region_profile,
//    fm_index_t* const fm_index,pattern_t* const pattern,bool* const allowed_enc,
//    const uint64_t region_length); // Much more todo
/*
 * Region Profile Generation (Fixed)
 */
GEM_INLINE void region_profile_generate_fixed(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,const uint64_t min_regions) {
  PROF_START(GP_REGION_PROFILE_ADAPTIVE);
  // Parameters
  const uint64_t rp_region_th = profile_model->region_th;
  const uint64_t rp_region_type_th = profile_model->region_type_th;
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  region_search_t* const regions = region_profile->filtering_region;
  const uint64_t max_region_length = key_length/min_regions;

  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    fprintf(stderr,"[%lu]",region_profile_num++);
    uint64_t i;
    for (i=0;i<key_length;++i) fprintf(stderr,"%c",dna_decode(key[i]));
    fprintf(stderr,"\n");
  }

  // Prepare rank_query
  const rank_mtable_t* const rank_mtable = fm_index->rank_table;
  rank_mquery_t rank_mquery;
  rank_mquery_new(&rank_mquery);

  // Initial region
  uint64_t num_regions = 0, num_standard_regions = 0, region_length = 0;
  uint64_t lo, hi, last_cut=0, max_regions=UINT64_MAX;
  uint64_t key_position = key_length;
  REGION_RESTART(key_position);
  while (key_position > 0) {
    --key_position;
    // Handling wildcards
    if (!allowed_enc[key[key_position]]) {
      while (key_position > 0 && !allowed_enc[key[key_position-1]]) --key_position;
      REGION_RESTART(key_position); region_length = 0; ++last_cut;
      continue;
    }

    // Rank query
    const uint8_t enc_char = key[key_position];
    if (!rank_mquery_is_exhausted(&rank_mquery)) {
      rank_mquery_add_char(rank_mtable,&rank_mquery,enc_char);
      rank_mtable_fetch(rank_mtable,&rank_mquery,&lo,&hi);
    } else {
      if (gem_expect_false(bwt_is_same_bucket(lo,hi))) {
        bwt_erank_interval(fm_index->bwt,enc_char,lo,hi,&lo,&hi);
      } else {
        lo = bwt_erank(fm_index->bwt,enc_char,lo);
        hi = bwt_erank(fm_index->bwt,enc_char,hi);
      }
    }
    ++region_length;

    // Determine region scope
    const uint64_t num_candidates = hi-lo;
    gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) { fprintf(stderr," %lu",num_candidates); }
    if (gem_expect_true(num_candidates <= rp_region_th || region_length >= max_region_length)) {
      REGION_CLOSE(key_position,lo,hi);
      REGION_RESTART(key_position); region_length = 0;
    }
  }
  // Check number of regions
  if (num_regions == 0) {
    if (regions[0].start == key_position) {
      regions[0].end = 0;
      regions[0].lo = lo;
      regions[0].hi = hi;
      region_profile->num_filtering_regions = 1;
      region_profile->num_standard_regions = 1;
    } else {
      region_profile->num_filtering_regions = 0;
      region_profile->num_standard_regions = 0;
    }
  } else {
    // We extend the last region
    region_profile->num_filtering_regions = num_regions;
    region_profile->num_standard_regions = num_standard_regions;
  }
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    fprintf(stderr,"\n");
  }
  PROF_STOP(GP_REGION_PROFILE_ADAPTIVE);
}
/*
 * Region Profile Adaptive
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
 *       - Regular Regions: Those regions which generates some few results in a exact search (> rp_region_type_th)
 *       - Unique Regions: Those regions which generates very few matches or none (<= rp_region_type_th)
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
GEM_INLINE void region_profile_generate_adaptive(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,
    const uint64_t max_regions,const bool allow_zero_regions) {
  PROF_START(GP_REGION_PROFILE_ADAPTIVE);
  const uint64_t rank_mquery_min_fetch_level = 0; // fm_index->proper_length/2; // FIXME => Change index
  const uint64_t rp_region_th = profile_model->region_th;
  const uint64_t rp_max_steps = profile_model->max_steps;
  const uint64_t rp_dec_factor = profile_model->dec_factor;
  const uint64_t rp_region_type_th = profile_model->region_type_th;
  region_search_t* const regions = region_profile->filtering_region;
  uint64_t num_regions = 0, num_standard_regions = 0, last_cut = 0;
  uint64_t lo, hi, hi_cut=0, lo_cut=0, expected_count, max_steps;

  // DEBUG
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    static uint64_t region_profile_num = 0;
    fprintf(stderr,"[%lu]",region_profile_num++);
    uint64_t i;
    for (i=0;i<key_length;++i) fprintf(stderr,"%c",dna_decode(key[i]));
    fprintf(stderr,"\n");
  }

  // Prepare rank_query
  const rank_mtable_t* const rank_mtable = fm_index->rank_table;
  rank_mquery_t rank_mquery;
  rank_mquery_new(&rank_mquery);

  // Initial region
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  uint64_t key_position = key_length;
  REGION_RESTART(key_position);
  while (key_position > 0) {
    --key_position;
    // Handling wildcards
    if (!allowed_enc[key[key_position]]) {
      if (last_cut != 0) REGION_CLOSE(key_position+1,lo,hi);
      while (key_position > 0 && !allowed_enc[key[key_position-1]]) --key_position;
      REGION_RESTART(key_position);
      continue;
    }

    // Rank query
    const uint8_t enc_char = key[key_position];
    if (!rank_mquery_is_exhausted(&rank_mquery)) {
      rank_mquery_add_char(rank_mtable,&rank_mquery,enc_char);
      if (rank_mquery.level > rank_mquery_min_fetch_level) {
        rank_mtable_fetch(rank_mtable,&rank_mquery,&lo,&hi);
      }
    } else {
      if (gem_expect_false(bwt_is_same_bucket(lo,hi))) {
        bwt_erank_interval(fm_index->bwt,enc_char,lo,hi,&lo,&hi);
      } else {
        lo = bwt_erank(fm_index->bwt,enc_char,lo);
        hi = bwt_erank(fm_index->bwt,enc_char,hi);
      }
    }
    // Determine region scope
    const uint64_t num_candidates = hi-lo;
    gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
      fprintf(stderr," %lu",num_candidates);
    }
    if (gem_expect_true(num_candidates > rp_region_th)) continue;

    // Check number of candidates (any or none)
    if (num_candidates > 0) {
      // We have candidate/s
      if (gem_expect_false(key_position == 0)) {
        // End-of-read
        REGION_CLOSE(key_position,lo,hi);
        REGION_RESTART(key_position);
      } else if (last_cut == 0) {
        // First Cut-Point
        REGION_SAVE_CUT_POINT(key_position,hi,lo);
        expected_count = num_candidates;
        max_steps = rp_max_steps;
      } else {
        // Check Cut-Point progress
        expected_count /= rp_dec_factor;
        if (num_candidates <= expected_count || num_candidates<=rp_region_type_th) {
          REGION_SAVE_CUT_POINT(key_position,hi,lo); // Refresh cut point
        }
        // Check maximum steps allowed for optimize region
        if (--max_steps == 0) {
          key_position = last_cut;
          REGION_CLOSE(key_position,lo_cut,hi_cut);
          REGION_RESTART(key_position);
        }
      }
    } else if (gem_expect_false(allow_zero_regions || last_cut == 0)) {
      // No candidates and/or no cut-point (we allow zero-regions)
      REGION_CLOSE(key_position,lo,hi);
      REGION_RESTART(key_position);
    } else {
      // We don't allow zero-regions (restore last cut-point)
      key_position = last_cut;
      REGION_CLOSE(key_position,lo_cut,hi_cut);
      REGION_RESTART(key_position);
    }
  }
  // Check number of regions
  if (num_regions == 0) {
    if (regions[0].start == key_length) {
      regions[0].end = 0;
      regions[0].lo = lo;
      regions[0].hi = hi;
      region_profile->num_filtering_regions = 1;
      region_profile->num_standard_regions = 1;
    } else {
      region_profile->num_filtering_regions = 0;
      region_profile->num_standard_regions = 0;
    }
  } else {
    // We extend the last region
    region_profile->num_filtering_regions = num_regions;
    region_profile->num_standard_regions = num_standard_regions;
    if (allow_zero_regions) {
      // NOTE Extending or not depending on the setup notably hits sensitivity
      region_profile_extend_last_region(region_profile,fm_index,key,allowed_enc,rp_region_type_th);
    }
  }
  gem_cond_debug_block(REGION_PROFILE_DEBUG_PRINT_PROFILE) {
    fprintf(stderr,"\n");
  }
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
 * Region Profile Utils
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
/*
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

*/




