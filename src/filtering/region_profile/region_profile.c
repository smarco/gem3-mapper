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
 *   Region-Profile data structure provides support to store a key partition
 *   for a candidate generation stage within a filtering process
 */

#include "align/pattern/pattern.h"
#include "filtering/region_profile/region_profile.h"
#include "fm_index/fm_index_search.h"

/*
 * Constants
 */
#define REGION_PROFILE_MIN_REGIONS_ALLOCATED 5

/*
 * Setup
 */
void region_profile_init(
    region_profile_t* const region_profile,
    const uint64_t pattern_length) {
  // Set length & max-expected regions
  region_profile->pattern_length = pattern_length;
  region_profile->max_expected_regions =
      MAX(REGION_PROFILE_MIN_REGIONS_ALLOCATED,DIV_CEIL(pattern_length,10));
  region_profile->filtering_region = NULL;
}
void region_profile_model_init(
    region_profile_model_t* const region_profile_model) {
  // Strategy
  region_profile_model->strategy = region_profile_adaptive;
  // General parameters (depending on the strategy)
  region_profile_model->num_regions = 5;
  region_profile_model->region_length = 20;
  region_profile_model->region_step = 20;
  region_profile_model->region_error = 0;
  region_profile_model->max_candidates = 100;
  // Adaptive parameters. LightWeight Scheme = (20,4,1)
  region_profile_model->region_th = 20;
  region_profile_model->max_steps = 4;
  region_profile_model->dec_factor = 1;
}
void region_profile_clear(
    region_profile_t* const region_profile) {
  region_profile->num_filtering_regions = 0;
  region_profile->num_filtered_regions = 0;
  region_profile->total_candidates = 0;
  region_profile->max_region_length = 0;
  region_profile->avg_region_length = 0;
  region_profile->kmer_frequency = 0.0;
}
void region_profile_inject_mm(
    region_profile_t* const region_profile,
    mm_allocator_t* const mm_allocator) {
  region_profile->mm_allocator = mm_allocator;
}
void region_profile_destroy(
    region_profile_t* const region_profile) {
  if (region_profile->filtering_region != NULL) {
    mm_allocator_free(region_profile->mm_allocator,region_profile->filtering_region);
  }
}
/*
 * Allocator
 */
void region_profile_allocate_regions(
    region_profile_t* const region_profile,
    const uint64_t num_regions) {
  if (num_regions > 0) {
    region_profile->filtering_region =
        mm_allocator_calloc(region_profile->mm_allocator,num_regions,region_search_t,false);
  } else {
    region_profile->filtering_region = NULL;
  }
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
 * Region Query
 */
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
      BWT_ERANK_TICK();
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
      BWT_ERANK_TICK();
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
      BWT_ERANK_TICK();
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
void region_profile_query_regions(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key) {
  // Iterate over all regions
  REGION_PROFILE_ITERATE(region_profile,region,i) {
    // Query
    fm_index_bsearch(fm_index,key+region->begin,
        region->end-region->begin,&region->hi,&region->lo);
  }
}
void region_profile_extend_last_region(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key) {
  // Tries to extend the last region of the profile to the end of the key (position 0)
  // Merge the tail with the last region
  region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
  if (last_region->begin != 0) {
    // Continue the search
    while (last_region->begin > 0) {
      if (last_region->hi==last_region->lo) break;
      // Query step
      const uint8_t enc_char = key[last_region->begin-1];
      if (enc_char == ENC_DNA_CHAR_N) break;
      last_region->lo=bwt_erank(fm_index->bwt,enc_char,last_region->lo);
      last_region->hi=bwt_erank(fm_index->bwt,enc_char,last_region->hi);
      --last_region->begin;
    }
    // Extend beyond zero interval // TODO
    // while (last_region->begin > 0 && allowed_enc[key[--last_region->begin]]);
  }
}
/*
 * Region Search Prepare
 */
void region_profile_fill_gaps_add(
    region_search_t* const filtering_region_filled,
    uint64_t* const num_filtering_regions_filled,
    const uint64_t begin_region,
    const uint64_t end_region,
    const uint64_t degree) {
  const uint64_t region_idx = (*num_filtering_regions_filled)++;
  filtering_region_filled[region_idx].begin = begin_region;
  filtering_region_filled[region_idx].end = end_region;
  filtering_region_filled[region_idx].degree = degree;
}
void region_profile_fill_gaps(
    region_profile_t* const region_profile,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t num_wildcards) {
  // Parameters
  mm_allocator_t* const mm_allocator = region_profile->mm_allocator;
  // Sort regions by position
  region_profile_sort_by_position(region_profile);
  // Current filtering-regions
  region_search_t* filtering_region = region_profile->filtering_region;
  uint64_t regions_left = region_profile->num_filtering_regions;
  // Allocate new filtering-regions
  const uint64_t max_regions = 2*(num_wildcards+regions_left)+1;
  region_search_t* const filtering_region_filled =
      mm_allocator_calloc(mm_allocator,MIN(key_length,max_regions),region_search_t,false);
  uint64_t num_filtering_regions_filled = 0;
  // Fill gaps
  uint64_t begin = 0, end = 0;
  for (;end<key_length;) {
    // Check current region profile
    if (regions_left>0 && filtering_region->begin==end) {
      // Fill gap
      if (begin < end) {
        region_profile_fill_gaps_add(filtering_region_filled,
            &num_filtering_regions_filled,begin,end,0);
      }
      region_profile_fill_gaps_add(filtering_region_filled,
          &num_filtering_regions_filled,filtering_region->begin,
          filtering_region->end,filtering_region->degree);
      // Restart
      begin = filtering_region->end;
      end   = filtering_region->end;
      // Next region
      --(regions_left);
      ++(filtering_region);
      continue;
    }
    // Next
    ++end;
  }
  // Fill gap
  if (begin < end) {
    region_profile_fill_gaps_add(filtering_region_filled,
        &num_filtering_regions_filled,begin,end,0);
  }
  // Set new region-profile
  region_profile->filtering_region = filtering_region_filled;
  region_profile->num_filtering_regions = num_filtering_regions_filled;
}
void region_profile_merge_small_regions(
    region_profile_t* const region_profile,
    const uint64_t proper_length) {
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  region_search_t* filtering_region = region_profile->filtering_region;
  region_search_t* filtering_region_before = NULL;
  region_search_t* filtering_region_after = NULL;
  region_search_t* filtering_region_out = region_profile->filtering_region;
  uint64_t i;
  for (i=0;i<num_filtering_regions;) {
    // Detect small regions bound not to precondition well
    const uint64_t region_length = filtering_region->end - filtering_region->begin;
    if (filtering_region->degree==0 && region_length < proper_length) {
      // Merge the region with the one before or after
      filtering_region_after = (i==num_filtering_regions-1) ? NULL : filtering_region+1;
      if (filtering_region_before && filtering_region_before->degree==0) {
        // Merge with region-before
        filtering_region_before->end = filtering_region->end;
        filtering_region_before->degree = 0;
        ++filtering_region; ++i;
        continue;
      } else if (filtering_region_after && filtering_region_after->degree==0 &&
                (filtering_region_after->end - filtering_region_after->begin) >= proper_length) {
        // Merge with region-after
        filtering_region_out->begin = filtering_region->begin;
        filtering_region_out->end = filtering_region_after->end;
        filtering_region_out->degree = 0;
        filtering_region_before = filtering_region_out;
        ++filtering_region_out;
        filtering_region += 2; i += 2;
        continue;
      } else if (filtering_region_before) {
        // Merge with region-before
        filtering_region_before->end = filtering_region->end;
        filtering_region_before->degree = 0;
        ++filtering_region; ++i;
        continue;
      } else if (filtering_region_after) {
        // Merge with region-after
        filtering_region_out->begin = filtering_region->begin;
        filtering_region_out->end = filtering_region_after->end;
        filtering_region_out->degree = 0;
        filtering_region_before = filtering_region_out;
        ++filtering_region_out;
        filtering_region += 2; i += 2;
        continue;
      }
    }
    // Copy region
    if (filtering_region_out != filtering_region) {
      *filtering_region_out = *filtering_region;
    }
    filtering_region_before = filtering_region_out;
    ++filtering_region_out;
    ++filtering_region; ++i;
  }
  // Set total number of regions
  region_profile->num_filtering_regions = filtering_region_out - region_profile->filtering_region;
}
uint64_t region_profile_compute_max_complete_strata(region_profile_t* const region_profile) {
  // Parameters
  const region_search_t* const filtering_region = region_profile->filtering_region;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  // Traverse all regions
  uint64_t i, max_accumulated_error = 0;
  for (i=0;i<num_filtering_regions;++i) {
    max_accumulated_error += filtering_region[i].degree;
  }
  return max_accumulated_error;
}
void region_profile_compute_error_limits(
    region_profile_t* const region_profile,
    const uint64_t max_complete_strata,
    const uint64_t max_search_error) {
  // Parameters
  region_search_t* const filtering_region = region_profile->filtering_region;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  // Traverse all regions
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    filtering_region[i].min = filtering_region[i].degree;
    filtering_region[i].max = max_search_error - (max_complete_strata-filtering_region[i].min);
  }
}
/*
 * kmer Frequency
 */
void region_profile_compute_kmer_frequency(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length) {
  // Parameters
  mm_allocator_t* const mm_allocator = region_profile->mm_allocator;
  // Push
  mm_allocator_push_state(mm_allocator);
  // Init mquery
  rank_mtable_t* const rank_mtable = fm_index->rank_table;
  const uint64_t max_samples = DIV_CEIL(key_length,RANK_MTABLE_SEARCH_DEPTH);
  rank_mquery_t* const rank_mquery = mm_allocator_calloc(mm_allocator,max_samples,rank_mquery_t,false);
  // Traverse the read & compute the frequency of contiguous kmers
  int64_t i, samples = 0;
  rank_mquery_t* current_rank_mquery = rank_mquery;
  rank_mquery_new(current_rank_mquery);
  for (i=key_length-1;i>=0;--i) {
    // Fetch character
    const uint8_t enc_char = key[i];
    if (enc_char == ENC_DNA_CHAR_N) {
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
  mm_allocator_pop_state(mm_allocator);
}
/*
 * Cmp
 */
int region_profile_cmp(
    region_profile_t* const region_profile_a,
    region_profile_t* const region_profile_b) {
  // Compare number of regions
  if (region_profile_a->num_filtering_regions !=
      region_profile_b->num_filtering_regions) return 1;
  // Compare regions
  const uint64_t num_regions = region_profile_a->num_filtering_regions;
  region_search_t* filtering_region_a = region_profile_a->filtering_region;
  region_search_t* filtering_region_b = region_profile_b->filtering_region;
  uint64_t i;
  for (i=0;i<num_regions;++i,++filtering_region_a,++filtering_region_b) {
    if (filtering_region_a->begin != filtering_region_b->begin) return 1;
    if (filtering_region_a->end != filtering_region_b->end) return 1;
    if (filtering_region_a->hi != filtering_region_b->hi) return 1;
    if (filtering_region_a->lo != filtering_region_b->lo) return 1;
  }
  // Return equals
  return 0;
}
/*
 * Sort
 */
int region_search_cmp_candidates(const region_search_t* const a,const region_search_t* const b) {
  return (int)(a->hi - a->lo) - (int)(b->hi - b->lo);
}
#define VECTOR_SORT_NAME                 region_search_by_candidates
#define VECTOR_SORT_TYPE                 region_search_t
#define VECTOR_SORT_CMP(a,b)             region_search_cmp_candidates(a,b)
#include "utils/vector_sort.h"
void region_profile_sort_by_candidates(region_profile_t* const region_profile) {
  buffer_sort_region_search_by_candidates(
      region_profile->filtering_region,region_profile->num_filtering_regions);
}
int region_search_cmp_position(const region_search_t* const a,const region_search_t* const b) {
  return (int)a->begin - (int)b->begin;
}
#define VECTOR_SORT_NAME                 region_search_by_position
#define VECTOR_SORT_TYPE                 region_search_t
#define VECTOR_SORT_CMP(a,b)             region_search_cmp_position(a,b)
#include "utils/vector_sort.h"
void region_profile_sort_by_position(region_profile_t* const region_profile) {
  buffer_sort_region_search_by_position(
      region_profile->filtering_region,region_profile->num_filtering_regions);
}
/*
 * Display
 */
void region_profile_print_region(
    FILE* const stream,
    region_search_t* const region,
    const uint64_t position,
    const bool display_error_limits) {
  if (!display_error_limits) {
    tab_fprintf(stream,
        "    [%"PRIu64"]\t"
        "region=[%"PRIu64",%"PRIu64")\t"
        "lo=%"PRIu64"\thi=%"PRIu64"\t"
        "cand=%"PRIu64"\n",
        position,region->begin,region->end,
        region->lo,region->hi,
        region->hi-region->lo);
  } else {
    tab_fprintf(stream,
        "    [%"PRIu64"]\t"
        "region=[%"PRIu64",%"PRIu64")\t"
        "error=[%"PRIu64",%"PRIu64"]\n",
        position,region->begin,region->end,
        region->min,region->max);
  }
}
void region_profile_print_region_pretty(
    FILE* const stream,
    const uint64_t region_idx,
    const uint64_t filter_degree,
    const uint64_t region_begin,
    const uint64_t region_end,
    const uint64_t region_candidates) {
  int64_t i;
  // Print candidates label
  const uint64_t region_length = region_end-region_begin;
  if (region_candidates!=UINT64_MAX) {
    tab_fprintf(stream,"                     ");
    const uint64_t ciphers = integer_num_ciphers(region_candidates);
    int64_t offset = region_begin + (region_length)/2 - ciphers/2;
    for (i=0;i<offset;++i) fprintf(stream," ");
    fprintf(stream,"%"PRIu64"\n",region_candidates);
  }
  // Print region
  tab_fprintf(stream,"    [#%03"PRIu64"]%s ",region_idx,filter_degree==0 ? "(Selected)" : "          ");
  for (i=0;i<region_begin;++i) fprintf(stream," ");
  if (region_length <= 1) {
    fprintf(stream,"|\n");
  } else {
    fprintf(stream,"|"); ++i;
    for (;i<region_end-1;++i) fprintf(stream,"-");
    fprintf(stream,"|\n");
  }
}
void region_profile_print(
    FILE* const stream,
    const region_profile_t* const region_profile,
    const bool display_error_limits) {
  tab_fprintf(stream,"[GEM]>Region.Profile\n");
  tab_fprintf(stream,"  => Pattern.length %"PRIu64"\n",region_profile->pattern_length);
  tab_fprintf(stream,"  => Num.Filtering.Regions %"PRIu64"\n",region_profile->num_filtering_regions);
  tab_fprintf(stream,"  => Filtering.Regions\n");
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    region_profile_print_region(stream,region,position,display_error_limits);
  }
  fflush(stream);
}
void region_profile_print_pretty(
    FILE* const stream,
    const region_profile_t* const region_profile,
    const char* const label,
    bool print_all_regions) {
  // Count selected regions
  tab_fprintf(stream,"[GEM]>Region.Profile (%s) (regions=%lu/%lu,cand=%lu)\n",
      label,region_profile->num_filtered_regions,
      region_profile->num_filtering_regions,
      region_profile->total_candidates);
  region_profile_print_region_pretty(stream,
      region_profile->num_filtering_regions,1,0,
      region_profile->pattern_length,UINT64_MAX);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    if (region->degree>0 || print_all_regions) {
      region_profile_print_region_pretty(
          stream,position,region->degree,
          region->begin,region->end,region->hi-region->lo);
    }
  }
  fflush(stream);
}

