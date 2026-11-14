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
static void region_profile_generator_setup(
    region_profile_generator_t* const generator,
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool allow_zero_regions,
    const bool allocate_regions) {
  // Region profile
  generator->region_profile = region_profile;
  region_profile_clear(region_profile);
  if (allocate_regions) {
    region_profile_allocate_regions(region_profile,region_profile->max_expected_regions); // Allocate
  }
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
void region_profile_generator_init(
    region_profile_generator_t* const generator,
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool allow_zero_regions) {
  region_profile_generator_setup(generator,
      region_profile,fm_index,key,key_length,allow_zero_regions,true);
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
 * Region Profile Adaptive computed ahead (several keys in lockstep)
 */
#define REGION_AHEAD_UNSET UINT64_MAX
static bool region_profile_generator_adaptive_step(
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model,
    const uint64_t proper_length) {
  // One iteration of region_profile_generate_adaptive's loop (false once it would have left it)
  if (generator->key_position == 0) return false;
  if (generator->region_profile->num_filtering_regions >=
      generator->region_profile->max_expected_regions) return false;
  --(generator->key_position);
  const uint8_t enc_char = generator->key[generator->key_position];
  if (enc_char == ENC_DNA_CHAR_N) {
    region_profile_generator_disallow_character(generator,profile_model);
  } else {
    region_profile_query_character(generator->fm_index,
        &generator->rank_mquery,&generator->lo,&generator->hi,enc_char);
    region_profile_generator_add_character(generator,profile_model,proper_length);
  }
  return true;
}
static void region_profile_generator_prefetch_step(
    const region_profile_generator_t* const generator,
    const rank_mtable_t* const rank_mtable,
    const bwt_t* const bwt) {
  // Prefetch what the next region_profile_query_character of this generator reads
  if (generator->key_position == 0) return;
  const uint8_t enc_char = generator->key[generator->key_position-1];
  if (enc_char == ENC_DNA_CHAR_N) return;
  if (!rank_mquery_is_exhausted(&generator->rank_mquery)) {
    const uint64_t level = generator->rank_mquery.level + 1;
    if (level >= rank_mtable->min_matching_depth) {
      const uint64_t hi_position = generator->rank_mquery.hi_position + enc_char*rank_mtable->level_skip[level];
      const uint64_t* const rank = rank_mtable->sa_ranks_levels[level] + (hi_position-1);
      PREFETCH(rank);
      PREFETCH(rank+1);
    }
  } else {
    const uint64_t lo = generator->lo, hi = generator->hi;
    const uint64_t* const block_hi = bwt->bwt_mem + (hi/64) * ((16 * 8 + 64 * 3 + 64) / 8 / 8);
    PREFETCH(block_hi);
    PREFETCH(block_hi+4);
    PREFETCH(bwt->mayor_counters + (hi / ((1 << 10) * 64)) * 8);
    if ((lo/64) != (hi/64)) {
      const uint64_t* const block_lo = bwt->bwt_mem + (lo/64) * ((16 * 8 + 64 * 3 + 64) / 8 / 8);
      PREFETCH(block_lo);
      PREFETCH(block_lo+4);
      PREFETCH(bwt->mayor_counters + (lo / ((1 << 10) * 64)) * 8);
    }
  }
}
bool region_profile_ahead_prepare(
    region_profile_ahead_t* const ahead,
    const uint64_t key_length) {
  // As region_profile_init; the regions are marked to tell the fields the generator writes
  ahead->ready = false;
  if (key_length == 0 || key_length > REGION_PROFILE_AHEAD_MAX_KEY_LENGTH) return false;
  region_profile_t* const region_profile = &ahead->region_profile;
  region_profile_init(region_profile,key_length);
  if (region_profile->max_expected_regions > REGION_PROFILE_AHEAD_MAX_REGIONS) return false;
  region_profile->filtering_region = ahead->regions;
  region_profile->mm_allocator = NULL;
  ahead->key_length = key_length;
  uint64_t i;
  for (i=0;i<=region_profile->max_expected_regions;++i) {
    region_search_t* const region = ahead->regions + i;
    region->begin = REGION_AHEAD_UNSET;
    region->end = REGION_AHEAD_UNSET;
    region->degree = REGION_AHEAD_UNSET;
    region->hi = REGION_AHEAD_UNSET;
    region->lo = REGION_AHEAD_UNSET;
  }
  return true;
}
void region_profile_ahead_generate(
    region_profile_ahead_t* const* const aheads,
    const uint64_t num_aheads,
    fm_index_t* const fm_index,
    const region_profile_model_t* const profile_model) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Parameters
  const uint64_t proper_length = fm_index_get_proper_length(fm_index);
  const rank_mtable_t* const rank_mtable = fm_index->rank_table;
  const bwt_t* const bwt = fm_index->bwt;
  // Init (the regions are already in place: nothing is allocated)
  region_profile_generator_t generators[REGION_PROFILE_AHEAD_MAX_KEYS];
  uint64_t active[REGION_PROFILE_AHEAD_MAX_KEYS];
  uint64_t i, num_active = 0;
  for (i=0;i<num_aheads;++i) {
    region_profile_ahead_t* const ahead = aheads[i];
    region_profile_generator_setup(generators+i,&ahead->region_profile,
        fm_index,ahead->key,ahead->key_length,false,false);
    region_profile_generator_prefetch_step(generators+i,rank_mtable,bwt);
    active[num_active++] = i;
  }
  // Delimit regions, one character of each key in turn (independent rank queries in flight)
  while (num_active > 0) {
    uint64_t j = 0;
    while (j < num_active) {
      region_profile_generator_t* const generator = generators + active[j];
      if (region_profile_generator_adaptive_step(generator,profile_model,proper_length)) {
        region_profile_generator_prefetch_step(generator,rank_mtable,bwt);
        ++j;
      } else {
        active[j] = active[--num_active];
      }
    }
  }
  // Close
  for (i=0;i<num_aheads;++i) {
    region_profile_generator_close_profile(generators+i,profile_model);
    // A write to the last slot lands past the allocation in region_profile_generate_adaptive: leave it to that
    region_profile_ahead_t* const ahead = aheads[i];
    const region_search_t* const past = ahead->regions + ahead->region_profile.max_expected_regions;
    ahead->fm_index = fm_index;
    ahead->region_th = profile_model->region_th;
    ahead->max_steps = profile_model->max_steps;
    ahead->dec_factor = profile_model->dec_factor;
    ahead->ready =
        past->begin == REGION_AHEAD_UNSET && past->end == REGION_AHEAD_UNSET &&
        past->degree == REGION_AHEAD_UNSET && past->hi == REGION_AHEAD_UNSET &&
        past->lo == REGION_AHEAD_UNSET;
  }
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
bool region_profile_ahead_restore(
    region_profile_ahead_t* const ahead,
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const region_profile_model_t* const profile_model) {
  // Replay only a profile of this very call: region_profile_generate_adaptive(...,ALL,false)
  if (!ahead->ready) return false;
  ahead->ready = false;
  if (ahead->fm_index != fm_index) return false;
  if (ahead->key_length != key_length) return false;
  if (region_profile->pattern_length != key_length) return false;
  if (region_profile->max_expected_regions != ahead->region_profile.max_expected_regions) return false;
  if (ahead->region_th != profile_model->region_th ||
      ahead->max_steps != profile_model->max_steps ||
      ahead->dec_factor != profile_model->dec_factor) return false;
  if (memcmp(ahead->key,key,key_length) != 0) return false;
  // As region_profile_generator_init (same allocation, in the same place)
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,region_profile->max_expected_regions);
  // The fields the generator wrote (the others keep what the allocation holds)
  uint64_t i;
  for (i=0;i<region_profile->max_expected_regions;++i) {
    const region_search_t* const src = ahead->regions + i;
    region_search_t* const dst = region_profile->filtering_region + i;
    if (src->begin != REGION_AHEAD_UNSET) dst->begin = src->begin;
    if (src->end != REGION_AHEAD_UNSET) dst->end = src->end;
    if (src->degree != REGION_AHEAD_UNSET) dst->degree = src->degree;
    if (src->hi != REGION_AHEAD_UNSET) dst->hi = src->hi;
    if (src->lo != REGION_AHEAD_UNSET) dst->lo = src->lo;
  }
  const region_profile_t* const profiled = &ahead->region_profile;
  region_profile->num_filtering_regions = profiled->num_filtering_regions;
  region_profile->num_filtered_regions = profiled->num_filtered_regions;
  region_profile->total_candidates = profiled->total_candidates;
  region_profile->max_region_length = profiled->max_region_length;
  region_profile->avg_region_length = profiled->avg_region_length;
  return true;
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
