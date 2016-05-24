/*
 * PROJECT: GEMMapper
 * FILE: region_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/region_profile.h"
#include "data_structures/pattern.h"

/*
 * Constants
 */
#define REGION_PROFILE_MIN_REGIONS_ALLOCATED 5

/*
 * Setup
 */
void region_profile_new(
    region_profile_t* const region_profile,
    const uint64_t pattern_length,
    mm_stack_t* const mm_stack) {
  // Allocate memory for the region profile
  region_profile->pattern_length = pattern_length;
  // Filtering regions
  region_profile->max_regions_allocated = MAX(REGION_PROFILE_MIN_REGIONS_ALLOCATED,DIV_CEIL(pattern_length,10));
  region_profile->filtering_region = mm_stack_calloc(mm_stack,region_profile->max_regions_allocated,region_search_t,false);
  // Locator for region sorting
  region_profile->loc = mm_stack_calloc(mm_stack,pattern_length,region_locator_t,false);
}
void region_profile_clear(region_profile_t* const region_profile) {
  // Reset
  region_profile->num_filtering_regions = 0;
  region_profile->errors_allowed = 0;
  region_profile->total_candidates = 0;
  region_profile->num_standard_regions = 0;
  region_profile->num_unique_regions = 0;
  region_profile->num_zero_regions = 0;
  region_profile->max_region_length = 0;
  region_profile->kmer_frequency = 0.0;
}
/*
 * Accessors
 */
uint64_t region_get_num_regions(region_profile_t* const region_profile) {
  // All regions (region_unique + region_standard + region_gap)
  return region_profile->num_filtering_regions;
}
bool region_profile_has_exact_matches(region_profile_t* const region_profile) {
  // Condition (sufficient but not necessary) to determine
  //   if there are exact matches (taken into account the first region)
  return (region_profile->num_filtering_regions==1) &&
      (region_profile->filtering_region[0].begin == 0 &&
       region_profile->filtering_region[0].end == region_profile->pattern_length &&
       (region_profile->filtering_region[0].hi-region_profile->filtering_region[0].lo) > 0);
}
/*
 * Utils
 */
void region_profile_compute_kmer_frequency(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    mm_stack_t* const mm_stack) {
  // Push
  mm_stack_push_state(mm_stack);
  // Init mquery
  rank_mtable_t* const rank_mtable = fm_index->rank_table;
  const uint64_t max_samples = DIV_CEIL(key_length,RANK_MTABLE_SEARCH_DEPTH);
  rank_mquery_t* const rank_mquery = mm_stack_calloc(mm_stack,max_samples,rank_mquery_t,false);
  // Traverse the read & compute the frequency of contiguous kmers
  int64_t i, samples = 0;
  rank_mquery_t* current_rank_mquery = rank_mquery;
  rank_mquery_new(current_rank_mquery);
  for (i=key_length-1;i>=0;--i) {
    // Fetch character
    const uint8_t enc_char = key[i];
    if (!allowed_enc[enc_char]) {
      rank_mquery_new(current_rank_mquery); // Reset
      continue;
    }
    // Query
    if (!rank_mquery_is_exhausted(current_rank_mquery)) {
      rank_mquery_add_char(rank_mtable,current_rank_mquery,enc_char);
    } else {
      rank_mtable_prefetch(rank_mtable,current_rank_mquery);
      ++samples; // Next
      ++current_rank_mquery;
      rank_mquery_new(current_rank_mquery); // Reset
    }
  }
  // Compute the kmer average frequency
  float frequency = 0.0;
  for (i=0;i<samples;++i) {
    // Account
    uint64_t hi, lo;
    rank_mtable_fetch(rank_mtable,rank_mquery+i,&lo,&hi);
    frequency += gem_loge((float)(hi-lo));
  }
  region_profile->kmer_frequency = (double)( frequency/((float)samples*gem_loge(4)) );
  // Free
  mm_stack_pop_state(mm_stack);
}
void region_profile_query_character(
    fm_index_t* const fm_index,
    rank_mquery_t* const rank_mquery,
    uint64_t* const lo,
    uint64_t* const hi,
    const uint8_t enc_char) {
  if (!rank_mquery_is_exhausted(rank_mquery)) {
    rank_mtable_t* const rank_mtable = fm_index->rank_table;
    const uint64_t min_matching_depth = rank_mtable->min_matching_depth;
    rank_mquery_add_char(rank_mtable,rank_mquery,enc_char);
    if (rank_mquery->level >= min_matching_depth) {
      rank_mtable_fetch(rank_mtable,rank_mquery,lo,hi);
    }
  } else {
    bwt_t* const bwt = fm_index->bwt;
    if (gem_expect_false((*lo/64) == (*hi/64) /*bwt_is_same_bucket(*lo,*hi)*/)) {
      /*
       * [MANUAL INLINE] bwt_erank_interval(bwt,enc_char,*lo,*hi,lo,hi); // Apologizes for doing this
       */
      const uint64_t block_pos = *hi / 64;
      const uint64_t block_mod = *hi % 64;
      const uint64_t* const mayor_counters = bwt->mayor_counters + (*hi / ((1 << 10) * 64)) * 8;
      const uint64_t* const block_mem = bwt->bwt_mem + block_pos * ((16 * 8 + 64 * 3 + 64) / 8 / 8);
      // Fetching Regular DNA Characters
      const uint64_t sum_counters = mayor_counters[enc_char] + ((uint16_t*) block_mem)[enc_char];
      const uint64_t bitmap = (block_mem[2] ^ xor_table_3[enc_char])
                            & (block_mem[3] ^ xor_table_2[enc_char])
                            & (block_mem[4] ^ xor_table_1[enc_char]);
      *hi = sum_counters + (POPCOUNT_64((bitmap & uint64_mask_ones[(block_mod)])));
      *lo = sum_counters + (POPCOUNT_64((bitmap & uint64_mask_ones[(*lo % 64)])));
    } else {
      /*
       * [MANUAL INLINE]  *lo = bwt_erank(bwt,enc_char,*lo); // Apologizes for doing this
       */
      const uint64_t block_pos_lo = *lo / 64;
      const uint64_t block_mod_lo = *lo % 64;
      const uint64_t* const mayor_counters_lo = bwt->mayor_counters + (*lo / ((1 << 10) * 64)) * 8;
      const uint64_t* const block_mem_lo = bwt->bwt_mem + block_pos_lo * ((16 * 8 + 64 * 3 + 64) / 8 / 8);
      // Calculate the exclusive rank for the given DNA character
      const uint64_t sum_counters_lo = mayor_counters_lo[enc_char] + ((uint16_t*) block_mem_lo)[enc_char];
      const uint64_t bitmap_lo = (block_mem_lo[2] ^ xor_table_3[enc_char])
                               & (block_mem_lo[3] ^ xor_table_2[enc_char])
                               & (block_mem_lo[4] ^ xor_table_1[enc_char]);
      // Return rank
      *lo = sum_counters_lo + POPCOUNT_64((bitmap_lo & uint64_mask_ones[(block_mod_lo)]));
      /*
       * [MANUAL INLINE]  *hi = bwt_erank(bwt,enc_char,*hi); // Apologizes for doing this
       */
      const uint64_t block_pos_hi = *hi / 64;
      const uint64_t block_mod_hi = *hi % 64;
      const uint64_t* const mayor_counters_hi = bwt->mayor_counters + (*hi / ((1 << 10) * 64)) * 8;
      const uint64_t* const block_mem_hi = bwt->bwt_mem + block_pos_hi * ((16 * 8 + 64 * 3 + 64) / 8 / 8);
      // Calculate the exclusive rank for the given DNA character
      const uint64_t sum_counters_hi = mayor_counters_hi[enc_char] + ((uint16_t*) block_mem_hi)[enc_char];
      const uint64_t bitmap_hi = (block_mem_hi[2] ^ xor_table_3[enc_char])
                               & (block_mem_hi[3] ^ xor_table_2[enc_char])
                               & (block_mem_hi[4] ^ xor_table_1[enc_char]);
      // Return rank
      *hi = sum_counters_hi + POPCOUNT_64((bitmap_hi & uint64_mask_ones[(block_mod_hi)]));
    }
  }
}
void region_profile_extend_last_region(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const bool* const allowed_enc,
    const uint64_t rp_region_type_th) {
  // Tries to extend the last region of the profile to the end of the key (position 0)
  // Merge the tail with the last region
  region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
  if (last_region->begin != 0) {
    // Continue the search
    while (last_region->begin > 0) {
      if (last_region->hi==last_region->lo) break;
      // Query step
      const uint8_t enc_char = key[last_region->begin-1];
      if (!allowed_enc[enc_char]) break;
      last_region->lo=bwt_erank(fm_index->bwt,enc_char,last_region->lo);
      last_region->hi=bwt_erank(fm_index->bwt,enc_char,last_region->hi);
      --last_region->begin;
    }
    // Extend beyond zero interval // TODO
    // while (last_region->begin > 0 && allowed_enc[key[--last_region->begin]]);
    // Adjust the region end (if needed) and type
    const uint64_t count = last_region->hi-last_region->lo;
    if (last_region->type==region_standard && count<=rp_region_type_th) {
      last_region->type = region_unique;
      --(region_profile->num_standard_regions);
    }
  }
}
/*
 * Sort
 */
int region_profile_locator_cmp(const region_locator_t* const a,const region_locator_t* const b) {
  return (int)a->value - (int)b->value;
}
//void region_profile_locator_sort(
//    region_locator_t* const loc,
//    const uint64_t num_regions) {
//  qsort(loc,num_regions,sizeof(region_locator_t),(int (*)(const void *,const void *))region_profile_locator_cmp);
//}
#define VECTOR_SORT_NAME                 region_profile_locator
#define VECTOR_SORT_TYPE                 region_locator_t
#define VECTOR_SORT_CMP(a,b)             region_profile_locator_cmp(a,b)
#include "utils/vector_sort.h"
void region_profile_locator_sort(region_locator_t* const loc,const uint64_t num_regions) {
  buffer_sort_region_profile_locator(loc,num_regions);
}
void region_profile_sort_by_estimated_mappability(region_profile_t* const region_profile) {
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
void region_profile_sort_by_candidates(region_profile_t* const region_profile) {
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
/*
 * Display
 */
void region_profile_print_region(
    FILE* const stream,
    region_search_t* const region,
    const uint64_t position) {
  tab_fprintf(stream,"    [%"PRIu64"]\ttype=%s\tregion=[%"PRIu64",%"PRIu64")\t"
      "lo=%"PRIu64"\thi=%"PRIu64"\tcand=%"PRIu64"\n",position,
      region->type==region_unique ? "region_unique" :
          (region->type==region_standard ? "region_standard" : "region_gap"),
      region->begin,region->end,region->lo,region->hi,region->hi-region->lo);
}
void region_profile_print(
    FILE* const stream,
    const region_profile_t* const region_profile,
    const bool sorted) {
  tab_fprintf(stream,"[GEM]>Region.Profile\n");
  tab_fprintf(stream,"  => Pattern.length %"PRIu64"\n",region_profile->pattern_length);
  tab_fprintf(stream,"  => Num.Filtering.Regions %"PRIu64"\n",region_profile->num_filtering_regions);
  tab_fprintf(stream,"  => Num.Standard.Regions %"PRIu64"\n",region_profile->num_standard_regions);
  tab_fprintf(stream,"  => Errors.allowed %"PRIu64"\n",region_profile->errors_allowed);
  tab_fprintf(stream,"  => Filtering.Regions\n");
  if (!sorted) {
    REGION_PROFILE_ITERATE(region_profile,region,position) {
      region_profile_print_region(stream,region,position);
    }
  } else {
    REGION_LOCATOR_ITERATE(region_profile,region,position) {
      region_profile_print_region(stream,region,position);
    }
  }
  fflush(stream);
}

