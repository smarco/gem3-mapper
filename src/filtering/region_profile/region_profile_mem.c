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
 *   Region-Profile module provides functions to generate an
 *   fixed/static key partition (with a fixed length size)
 */

#include "filtering/region_profile/region_profile_mem.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "fm_index/fm_index_query.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * MEMs region profile
 */
void region_profile_compute_mem(
    region_search_t* const filtering_region,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t end_position) {
  // Init
  uint64_t prev_lo=0, lo = 0;
  uint64_t prev_hi=0, hi = fm_index_get_length(fm_index);
  int64_t position = end_position-1;
  while (position>=0) {
    // Fetch character
    const uint8_t enc_char = key[position];
    if (enc_char == ENC_DNA_CHAR_N) break;
    // Query index
    lo = bwt_sbasic_erank(fm_index->bwt,enc_char,lo);
    hi = bwt_sbasic_erank(fm_index->bwt,enc_char,hi);
    if (hi==lo) break;
    // Store
    prev_lo = lo;
    prev_hi = hi;
    // Next
    --position;
  }
  // Configure region
  filtering_region->begin = position+1;
  filtering_region->end = end_position;
  filtering_region->hi = prev_hi;
  filtering_region->lo = prev_lo;
}
void region_profile_generate_mem(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_candidates) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Init
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,key_length); // Allocate
  // Compute MEM for each key position
  region_search_t* filtering_region = region_profile->filtering_region; // Filtering regions
  uint64_t num_regions = 0;
  uint64_t i;
  for (i=key_length;i>0;) {
    // Fetch character
    const uint8_t enc_char = key[i-1];
    if (enc_char != ENC_DNA_CHAR_N) {
      // Compute MEM
      region_profile_compute_mem(filtering_region,fm_index,key,key_length,i);
      // This cond. is commented as otherwise it will generate SMEMs only
      // if (filtering_region->begin < last_begin_position) {
      // Next region
      ++filtering_region;
      ++num_regions;
    }
    // Next position
    --i;
  }
  region_profile->num_filtering_regions = num_regions;
  // Schedule to filter all
  region_profile_schedule_exact(region_profile,max_candidates);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
/*
 * SMEMs region profile
 */
typedef struct {
  uint64_t begin;
  uint64_t end;
  fm_2interval_t fm_2interval;
} smem_interval_t;
typedef struct {
  smem_interval_t* smem_interval;
  uint64_t num_smem_intervals;
} smem_interval_vector_t;
void smem_interval_vector_append(
    smem_interval_vector_t* const smem_interval_vector,
    fm_2interval_t* const interval,
    const uint64_t begin_position,
    const uint64_t end_position) {
  smem_interval_t* const smem_interval =
      smem_interval_vector->smem_interval + smem_interval_vector->num_smem_intervals;
  smem_interval->fm_2interval = *interval;
  smem_interval->begin = begin_position;
  smem_interval->end = end_position;
  ++(smem_interval_vector->num_smem_intervals);
}
uint64_t region_profile_compute_smem(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t init_position,
    smem_interval_vector_t* const smem_interval_vector) {
  // Init
  smem_interval_vector->num_smem_intervals = 0;
  fm_2interval_t curr_interval, next_interval;
  fm_index_2query_init(fm_index,&curr_interval);
  // Right extend
  int64_t i, idx, next_position = init_position+1;
  for (i=init_position;i<key_length;++i) {
    // Fetch character
    const uint8_t enc_char = key[i];
    if (enc_char == ENC_DNA_CHAR_N) {
      next_interval.backward_hi = 0;
      next_interval.backward_lo = 0;
    } else {
      // Extend forward
      fm_index_2query_forward_query(fm_index,&curr_interval,&next_interval,enc_char);
    }
    // Check occurrences
    const uint64_t next_occ = next_interval.backward_hi-next_interval.backward_lo;
    if (next_occ==0) break;
    const uint64_t curr_occ = curr_interval.backward_hi-curr_interval.backward_lo;
    if (i > init_position && curr_occ != next_occ) {
      smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
    }
    // Next
    next_position = i+1;
    SWAP(curr_interval,next_interval); // Swap
  }
  smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
  // Left extend
  const uint64_t num_smem_intervals = smem_interval_vector->num_smem_intervals;
  uint64_t leftmost = init_position+1;
  for (idx=num_smem_intervals-1;idx>=0;--idx) {
    smem_interval_t* const smem_interval = smem_interval_vector->smem_interval + idx;
    for (i=init_position-1;i>=0;--i) {
      // Fetch character
      const uint8_t enc_char = key[i];
      if (enc_char == ENC_DNA_CHAR_N) {
        next_interval.backward_hi = 0;
        next_interval.backward_lo = 0;
        break;
      } else {
        // Extend forward
        fm_index_2query_backward_query(fm_index,
            &smem_interval->fm_2interval,&next_interval,enc_char);
      }
      // Check occurrences
      const uint64_t next_occ = next_interval.backward_hi - next_interval.backward_lo;
      if (next_occ==0) break;
      // Update
      smem_interval->begin = i;
      smem_interval->fm_2interval = next_interval;
    }
    // Check leftmost position
    if (smem_interval->begin < leftmost) {
      // Configure region
      region_search_t* const filtering_region =
          region_profile->filtering_region + region_profile->num_filtering_regions;
      filtering_region->begin = smem_interval->begin;
      filtering_region->end = smem_interval->end;
      filtering_region->hi = smem_interval->fm_2interval.backward_hi;
      filtering_region->lo = smem_interval->fm_2interval.backward_lo;
      ++(region_profile->num_filtering_regions);
      // Update
      leftmost = smem_interval->begin;
    }
  }
  // Return
  return next_position;
}
void region_profile_generate_smem(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_candidates) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Parameters
  mm_allocator_t* const mm_allocator = region_profile->mm_allocator;
  // Init
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,key_length);
  // Allocate
  mm_allocator_push_state(mm_allocator);
  smem_interval_vector_t smem_interval_vector;
  smem_interval_vector.smem_interval =
      mm_allocator_calloc(mm_allocator,key_length,smem_interval_t,false);
  // Compute SMEMs for each key position
  uint64_t i;
  region_profile->num_filtering_regions = 0;
  for (i=0;i<key_length;) {
    // Compute SMEMs
    if (key[i] == ENC_DNA_CHAR_N) {
      ++i;
    } else {
      const uint64_t next_pos = region_profile_compute_smem(
          region_profile,fm_index,key,key_length,i,&smem_interval_vector);
      i = next_pos;
    }
  }
  // Schedule to filter all
  region_profile_schedule_exact(region_profile,max_candidates);
  // Free
  mm_allocator_pop_state(mm_allocator);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
///*
// * PAPER IMPLEMENTATION (correcting a couple of errors)
// */
//void region_profile_compute_smem(
//    region_profile_t* const region_profile,
//    fm_index_t* const fm_index,
//    const uint8_t* const key,
//    const uint64_t key_length,
//    const bool* const allowed_enc,
//    const uint64_t init_position,
//    smem_interval_vector_t* const smem_interval_vector) {
//  // Init
//  smem_interval_vector->num_smem_intervals = 0;
//  fm_2interval_t curr_interval, next_interval;
//  fm_index_2query_init(fm_index,&curr_interval);
//  // Right extend
//  int64_t i;
//  for (i=init_position;i<key_length;++i) {
//    // Fetch character
//    const uint8_t enc_char = key[i];
//    if (!allowed_enc[enc_char]) {
//      next_interval.backward_hi = 0;
//      next_interval.backward_lo = 0;
//    } else {
//      // Extend forward
//      fm_index_2query_forward_query(fm_index,&curr_interval,&next_interval,enc_char);
//    }
//    // Check occurrences
//    const uint64_t next_occ = next_interval.backward_hi-next_interval.backward_lo;
//    if (next_occ==0) break;
//    const uint64_t curr_occ = curr_interval.backward_hi-curr_interval.backward_lo;
//    if (i > init_position && curr_occ != next_occ) {
//      smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
//    }
//    // Swap
//    SWAP(curr_interval,next_interval);
//  }
//  smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
//  // Left extend
//  uint64_t idx, lowest_pos = key_length;
//  for (i=init_position-1;i>=0;--i) {
//    uint64_t lowest_occ = UINT64_MAX;
//    const uint64_t num_smem_intervals = smem_interval_vector->num_smem_intervals;
//    smem_interval_t* smem_interval_in = smem_interval_vector->smem_interval;
//    uint64_t num_intervals_out = 0;
//    for (idx=0;idx<num_smem_intervals;++idx,++smem_interval_in) {
//      // Fetch character
//      const uint8_t enc_char = key[i];
//      if (!allowed_enc[enc_char]) {
//        next_interval.backward_hi = 0;
//        next_interval.backward_lo = 0;
//      } else {
//        // Extend forward
//        fm_index_2query_backward_query(fm_index,
//            &smem_interval_in->fm_2interval,&next_interval,enc_char);
//      }
//      // Check occurrences
//      const uint64_t next_occ = next_interval.backward_hi - next_interval.backward_lo;
//      if (next_occ==0) {
//        if (num_intervals_out==0 && i+1<lowest_pos+1) {
//          lowest_pos = i;
//          // Add region current
//          region_search_t* const filtering_region =
//              region_profile->filtering_region + region_profile->num_filtering_regions;
//          filtering_region->begin = smem_interval_in->begin;
//          filtering_region->end = smem_interval_in->end;
//          filtering_region->hi = smem_interval_in->fm_2interval.backward_hi;
//          filtering_region->lo = smem_interval_in->fm_2interval.backward_lo;
//          ++(region_profile->num_filtering_regions);
//        }
//      } else {
//        if (lowest_occ!=next_occ) {
//          lowest_occ = next_occ;
//          smem_interval_vector->smem_interval[num_intervals_out].end = smem_interval_in->end;
//          smem_interval_vector->smem_interval[num_intervals_out].begin = i;
//          smem_interval_vector->smem_interval[num_intervals_out].fm_2interval = next_interval;
//          ++num_intervals_out;
//        }
//      }
//    }
//    // Set vector used
//    smem_interval_vector->num_smem_intervals = num_intervals_out;
//    // Check vector
//    if (num_intervals_out==0) break;
//  }
//  // Add remaining intervals
//  const uint64_t num_smem_intervals = smem_interval_vector->num_smem_intervals;
//  uint64_t num_intervals_out = 0;
//  for (idx=0;idx<num_smem_intervals;++idx) {
//    smem_interval_t* const smem_interval_in = smem_interval_vector->smem_interval + idx;
//    // Check occurrences
//    if (num_intervals_out==0) {
//      // Add region current
//      region_search_t* const filtering_region =
//          region_profile->filtering_region + region_profile->num_filtering_regions;
//      filtering_region->begin = smem_interval_in->begin;
//      filtering_region->end = smem_interval_in->end;
//      filtering_region->hi = smem_interval_in->fm_2interval.backward_hi;
//      filtering_region->lo = smem_interval_in->fm_2interval.backward_lo;
//      ++(region_profile->num_filtering_regions);
//      break;
//    }
//  }
//}

///*
// * PAPER FIXED
// */
//uint64_t region_profile_compute_smem(
//    region_profile_t* const region_profile,
//    fm_index_t* const fm_index,
//    const uint8_t* const key,
//    const uint64_t key_length,
//    const bool* const allowed_enc,
//    const uint64_t init_position,
//    smem_interval_vector_t* const smem_interval_vector) {
//  // Init
//  smem_interval_vector->num_smem_intervals = 0;
//  fm_2interval_t curr_interval, next_interval;
//  fm_index_2query_init(fm_index,&curr_interval);
//  // Right extend
//  int64_t i, next_position = init_position+1;
//  for (i=init_position;i<key_length;++i) {
//    // Fetch character
//    const uint8_t enc_char = key[i];
//    if (!allowed_enc[enc_char]) {
//      next_interval.backward_hi = 0;
//      next_interval.backward_lo = 0;
//    } else {
//      // Extend forward
//      fm_index_2query_forward_query(fm_index,&curr_interval,&next_interval,enc_char);
//    }
//    // Check occurrences
//    const uint64_t next_occ = next_interval.backward_hi-next_interval.backward_lo;
//    if (next_occ==0) break;
//    const uint64_t curr_occ = curr_interval.backward_hi-curr_interval.backward_lo;
//    if (i > init_position && curr_occ != next_occ) {
//      smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
//    }
//    // Next
//    next_position = i+1;
//    SWAP(curr_interval,next_interval); // Swap
//  }
//  smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
//  // Left extend
//  uint64_t idx;
//  for (i=init_position-1;i>=0;--i) {
//    uint64_t lowest_occ = UINT64_MAX;
//    const uint64_t num_smem_intervals = smem_interval_vector->num_smem_intervals;
//    smem_interval_t* smem_interval_in = smem_interval_vector->smem_interval;
//    uint64_t num_intervals_out = 0;
//    for (idx=0;idx<num_smem_intervals;++idx,++smem_interval_in) {
//      // Fetch character
//      const uint8_t enc_char = key[i];
//      if (!allowed_enc[enc_char]) {
//        next_interval.backward_hi = 0;
//        next_interval.backward_lo = 0;
//      } else {
//        // Extend forward
//        fm_index_2query_backward_query(fm_index,
//            &smem_interval_in->fm_2interval,&next_interval,enc_char);
//      }
//      // Check occurrences
//      const uint64_t next_occ = next_interval.backward_hi - next_interval.backward_lo;
//      if (next_occ==0) {
//        // Add last smem-interval
//        smem_interval_t* const smem_interval_last =
//            smem_interval_vector->smem_interval + (num_smem_intervals-1);
//        // Configure region
//        region_search_t* const filtering_region =
//            region_profile->filtering_region + region_profile->num_filtering_regions;
//        filtering_region->begin = smem_interval_last->begin;
//        filtering_region->end = smem_interval_last->end;
//        filtering_region->hi = smem_interval_last->fm_2interval.backward_hi;
//        filtering_region->lo = smem_interval_last->fm_2interval.backward_lo;
//        ++(region_profile->num_filtering_regions);
//        break;
//      } else {
//        if (lowest_occ==next_occ && num_intervals_out>0) --num_intervals_out;
//        lowest_occ = next_occ;
//        smem_interval_vector->smem_interval[num_intervals_out].end = smem_interval_in->end;
//        smem_interval_vector->smem_interval[num_intervals_out].begin = i;
//        smem_interval_vector->smem_interval[num_intervals_out].fm_2interval = next_interval;
//        ++num_intervals_out;
//      }
//    }
//    // Set vector used
//    smem_interval_vector->num_smem_intervals = num_intervals_out;
//    if (num_intervals_out==0) break;
//  }
//  // Add last interval
//  if (smem_interval_vector->num_smem_intervals > 0) {
//    // Add last smem-interval
//    smem_interval_t* const smem_interval_last =
//        smem_interval_vector->smem_interval + (smem_interval_vector->num_smem_intervals-1);
//    // Configure region
//    region_search_t* const filtering_region =
//        region_profile->filtering_region + region_profile->num_filtering_regions;
//    filtering_region->begin = smem_interval_last->begin;
//    filtering_region->end = smem_interval_last->end;
//    filtering_region->hi = smem_interval_last->fm_2interval.backward_hi;
//    filtering_region->lo = smem_interval_last->fm_2interval.backward_lo;
//    ++(region_profile->num_filtering_regions);
//  }
//  // Return
//  return next_position;
//}

///*
// * SMEM Implementation keeping MAX-MEM in case of conflict
// */
//uint64_t region_profile_compute_smem(
//    region_profile_t* const region_profile,
//    fm_index_t* const fm_index,
//    const uint8_t* const key,
//    const uint64_t key_length,
//    const bool* const allowed_enc,
//    const uint64_t init_position,
//    smem_interval_vector_t* const smem_interval_vector) {
//  // Init
//  smem_interval_vector->num_smem_intervals = 0;
//  fm_2interval_t curr_interval, next_interval;
//  fm_index_2query_init(fm_index,&curr_interval);
//  // Right extend
//  int64_t i;
//  for (i=init_position;i<key_length;++i) {
//    // Fetch character
//    const uint8_t enc_char = key[i];
//    if (!allowed_enc[enc_char]) {
//      next_interval.backward_hi = 0;
//      next_interval.backward_lo = 0;
//    } else {
//      // Extend forward
//      fm_index_2query_forward_query(fm_index,&curr_interval,&next_interval,enc_char);
//    }
//    // Check occurrences
//    const uint64_t next_occ = next_interval.backward_hi-next_interval.backward_lo;
//    if (next_occ==0) break;
//    const uint64_t curr_occ = curr_interval.backward_hi-curr_interval.backward_lo;
//    if (i > init_position && curr_occ != next_occ) {
//      smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
//    }
//    // Swap
//    SWAP(curr_interval,next_interval);
//  }
//  smem_interval_vector_append(smem_interval_vector,&curr_interval,init_position,i);
//  // Left extend
//  region_search_t largest_filtering_region;
//  largest_filtering_region.begin = 0;
//  largest_filtering_region.end = 0;
//  uint64_t idx;
//  for (i=init_position-1;i>=0;--i) {
//    uint64_t lowest_occ = UINT64_MAX;
//    const uint64_t num_smem_intervals = smem_interval_vector->num_smem_intervals;
//    smem_interval_t* smem_interval_in = smem_interval_vector->smem_interval;
//    uint64_t num_intervals_out = 0;
//    for (idx=0;idx<num_smem_intervals;++idx,++smem_interval_in) {
//      // Fetch character
//      const uint8_t enc_char = key[i];
//      if (!allowed_enc[enc_char]) {
//        next_interval.backward_hi = 0;
//        next_interval.backward_lo = 0;
//      } else {
//        // Extend forward
//        fm_index_2query_backward_query(fm_index,
//            &smem_interval_in->fm_2interval,&next_interval,enc_char);
//      }
//      // Check occurrences
//      const uint64_t next_occ = next_interval.backward_hi - next_interval.backward_lo;
//      if (next_occ==0) {
//        // Add last smem-interval (largest zero)
//        smem_interval_t* const smem_interval_last =
//            smem_interval_vector->smem_interval + (num_smem_intervals-1);
//        if (largest_filtering_region.end-largest_filtering_region.begin <
//            smem_interval_last->end-smem_interval_last->begin) {
//          largest_filtering_region.end = smem_interval_last->end;
//          largest_filtering_region.begin = smem_interval_last->begin;
//          largest_filtering_region.hi = smem_interval_last->fm_2interval.backward_hi;
//          largest_filtering_region.lo = smem_interval_last->fm_2interval.backward_lo;
//        }
//        break;
//      } else {
//        if (lowest_occ==next_occ && num_intervals_out>0) --num_intervals_out;
//        lowest_occ = next_occ;
//        smem_interval_vector->smem_interval[num_intervals_out].end = smem_interval_in->end;
//        smem_interval_vector->smem_interval[num_intervals_out].begin = i;
//        smem_interval_vector->smem_interval[num_intervals_out].fm_2interval = next_interval;
//        ++num_intervals_out;
//      }
//    }
//    // Set vector used
//    smem_interval_vector->num_smem_intervals = num_intervals_out;
//    if (num_intervals_out==0) break;
//  }
//  // Add last smem-interval (largest zero)
//  const uint64_t num_smem_intervals = smem_interval_vector->num_smem_intervals;
//  if (num_smem_intervals > 0) {
//    smem_interval_t* const smem_interval_last =
//        smem_interval_vector->smem_interval + (num_smem_intervals-1);
//    if (largest_filtering_region.end-largest_filtering_region.begin < smem_interval_last->end) {
//      largest_filtering_region.begin = smem_interval_last->begin;
//      largest_filtering_region.end = smem_interval_last->end;
//      largest_filtering_region.hi = smem_interval_last->fm_2interval.backward_hi;
//      largest_filtering_region.lo = smem_interval_last->fm_2interval.backward_lo;
//    }
//  }
//  // Configure region
//  region_profile->filtering_region[region_profile->num_filtering_regions] = largest_filtering_region;
//  ++(region_profile->num_filtering_regions);
//  // Return
//  return largest_filtering_region.end;
//}

