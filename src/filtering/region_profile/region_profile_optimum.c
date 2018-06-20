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
 */

#include "filtering/region_profile/region_profile_optimum.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "align/pattern/pattern.h"
#include "fm_index/fm_index_search.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Data structures
 */
typedef struct {
  uint64_t hi;
  uint64_t lo;
} optimum_table_interval_t;
typedef struct {
  bool region_taken;
  uint64_t candidates;
} optimum_fixed_cell_t;
typedef struct {
  uint64_t occ;
  uint64_t split_point;
} optimum_variable_cell_t;

/*
 * Optimal prefix selection: select those fixed-length
 *   regions with fewer candidates from all possible
 *   combinations (Dynamic region profile)
 */
uint64_t region_profile_generate_ops_fetch_region(
    optimum_table_interval_t* const optimum_table_interval,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t position,
    const uint64_t region_length) {
  // Check boundaries
  if (position+1 < region_length) {
    return UINT64_MAX;
  } else {
    // Check computed
    if (optimum_table_interval[position].hi==UINT64_MAX) {
      const uint64_t begin_pos = position-region_length+1;
      const uint64_t end_pos = position+1;
      fm_index_bsearch(fm_index,key+begin_pos,end_pos-begin_pos,
          &optimum_table_interval[position].hi,&optimum_table_interval[position].lo);
    }
    // Return number of candidates
    return optimum_table_interval[position].hi - optimum_table_interval[position].lo;
  }
}
void region_profile_generate_optimum_fixed(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t region_length,
    const uint64_t num_regions,
    const uint64_t max_candidates) {
  // Base check
  if (num_regions==0 || key_length/num_regions < region_length) {
    gem_warn_msg("OPS. Not enough bases to extract %"PRIu64" regions of %"PRIu64" bases",num_regions,region_length);
    region_profile->total_candidates = 0;
    region_profile->num_filtering_regions = 0;
    return;
  }
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Parameters
  mm_allocator_t* const mm_allocator = region_profile->mm_allocator;
  // Init
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,num_regions); // Allocate
  // Allocate
  mm_allocator_push_state(mm_allocator);
  optimum_table_interval_t* const optimum_table_interval =
      mm_allocator_calloc(mm_allocator,key_length,optimum_table_interval_t,false);
  optimum_fixed_cell_t* const optimum_fixed_cell =
      mm_allocator_calloc(mm_allocator,num_regions*key_length,optimum_fixed_cell_t,false);
  // Init
  int64_t i, level, current_idx;
  for (i=0;i<key_length;++i) {
    optimum_table_interval[i].hi = UINT64_MAX;
  }
  // Compute optimization vector
  for (level=0;level<num_regions;++level) {
    for (i=0;i<key_length;++i) {
      current_idx = level*key_length + i;
      // Compute cost of skip
      const uint64_t cost_skip = (i>0) ? optimum_fixed_cell[current_idx-1].candidates : UINT64_MAX;
      // Compute cost of region
      const uint64_t cost_candidates =
          region_profile_generate_ops_fetch_region(
              optimum_table_interval,fm_index,key,i,region_length);
      // Compute cost of taken
      uint64_t cost_taken;
      if (cost_candidates==UINT64_MAX || level==0) {
        cost_taken = cost_candidates;
      } else {
        if (i>=region_length) {
          const uint64_t prev_idx = (level-1)*key_length + i-region_length;
          if (optimum_fixed_cell[prev_idx].candidates == UINT64_MAX) {
            cost_taken = UINT64_MAX;
          } else {
            cost_taken = optimum_fixed_cell[prev_idx].candidates + cost_candidates;
          }
        } else {
          cost_taken = UINT64_MAX;
        }
      }
      // Compute minimum
      if (cost_skip <= cost_taken) {
        optimum_fixed_cell[current_idx].candidates = cost_skip;
        optimum_fixed_cell[current_idx].region_taken = false;
      } else {
        optimum_fixed_cell[current_idx].candidates = cost_taken;
        optimum_fixed_cell[current_idx].region_taken = true;
      }
    }
  }
  // Backtrace result from optimum: optimum_fixed_cell[num_regions-1][key_length-1]
  current_idx = (num_regions-1)*key_length + key_length-1;
  if (optimum_fixed_cell[current_idx].candidates == UINT64_MAX) {
    gem_warn_msg("OPS couldn't converge (not enough bases)");
    region_profile->total_candidates = 0;
    region_profile->num_filtering_regions = 0;
  } else {
    region_search_t* filtering_region = region_profile->filtering_region;
    region_profile->num_filtering_regions = 0;
    level = num_regions-1;
    i = key_length-1;
    while (level >= 0) {
      current_idx = level*key_length + i;
      if (!optimum_fixed_cell[current_idx].region_taken) {
        --i; // Skip
      } else {
        filtering_region->begin = i-region_length+1;
        filtering_region->end = i+1;
        filtering_region->hi = optimum_table_interval[i].hi;
        filtering_region->lo = optimum_table_interval[i].lo;
        ++(filtering_region);
        ++(region_profile->num_filtering_regions);
        // Next
        --level;
        i = i-region_length;
      }
    }
  }
  // Schedule to filter all
  region_profile_schedule_exact(region_profile,max_candidates);
  // Free
  mm_allocator_pop_state(mm_allocator);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
/*
 * Region Profile Optimum (Variable region-length)
 */
void region_profile_optimum_compute_occ(
    uint64_t** const occ,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length) {
  int64_t end, begin;
  // Starting from e, search backwards
  for (end=key_length-1;end>=0;--end) {
    // Fill column occ(,end)
    uint64_t lo = 0;
    uint64_t hi = fm_index_get_length(fm_index);
    for (begin=end;begin>=0;--begin) {
      const uint8_t enc_char = key[begin];
      if (enc_char == ENC_DNA_CHAR_N) {
        for (;begin>=0;--begin) occ[begin][end] = 0;
        break;
      }
      lo = bwt_sbasic_erank(fm_index->bwt,enc_char,lo);
      hi = bwt_sbasic_erank(fm_index->bwt,enc_char,hi);
      occ[begin][end] = hi-lo;
      if (occ[begin][end]==0) {
        for (;begin>=0;--begin) occ[begin][end] = 0;
        break;
      }
    }
  }
}
void region_profile_generate_optimum_variable(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t num_regions) {
  // Base check
  if (num_regions==0 || num_regions > key_length) {
    gem_warn_msg("Not enough bases to extract %"PRIu64" regions",num_regions);
    region_profile->total_candidates = 0;
    region_profile->num_filtering_regions = 0;
    return;
  }
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Parameters
  mm_allocator_t* const mm_allocator = region_profile->mm_allocator;
  // Init
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,num_regions); // Allocate
  // Compute occurrences table
  uint64_t i;
  uint64_t** const occ = mm_allocator_calloc(mm_allocator,key_length,uint64_t*,false);
  for (i=0;i<key_length;++i) {
    occ[i] = mm_allocator_calloc(mm_allocator,key_length,uint64_t,true);
  }
  region_profile_optimum_compute_occ(occ,fm_index,key,key_length);
  // Allocate optimum pattern partition (opp)
  const uint64_t last_position = key_length-1;
  optimum_variable_cell_t** const opp =
      mm_allocator_calloc(mm_allocator,num_regions,optimum_variable_cell_t*,true);
  for (i=0;i<num_regions;++i) {
    opp[i] = mm_allocator_calloc(mm_allocator,key_length,optimum_variable_cell_t,true);
  }
  // Init
  for (i=0;i<key_length;++i) {
    opp[0][i].occ = occ[i][last_position];
    opp[0][i].split_point = key_length;
  }
  // Compute the optimum
  uint64_t j, p, split_point;
  for (i=1;i<num_regions;++i) {
    // Find the minimum candidates split-point from each position
    for (j=0;j<key_length;++j) {
      split_point = 0;
      uint64_t split_point_occ = UINT32_MAX;
      for (p=j+1;p<key_length;++p) {
        const uint64_t total_occ = occ[j][p-1] + opp[i-1][p].occ;
        if (total_occ < split_point_occ) {
          split_point = p;
          split_point_occ = total_occ;
        }
      }
      opp[i][j].split_point = split_point;
      opp[i][j].occ = split_point_occ;
    }
  }
  // Compose optimum region profile
  region_profile->total_candidates = opp[num_regions-1][0].occ; // Minimum
  region_profile->num_filtering_regions = num_regions;
  p=0; split_point = 0;
  for (i=num_regions;i>0;--i) {
    split_point = opp[i-1][p].split_point;
    region_search_t* const filtering_region = &region_profile->filtering_region[num_regions-i];
    filtering_region->begin = p;
    filtering_region->end = split_point;
    filtering_region->lo = 0; //  FIXME  FIXME  FIXME  FIXME  FIXME  FIXME
    filtering_region->hi = occ[p][split_point-1];
    p = split_point;
  }
  // Schedule to filter all
  region_profile_schedule_exact(region_profile,ALL);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}

