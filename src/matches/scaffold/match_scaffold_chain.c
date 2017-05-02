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

#include "matches/scaffold/match_scaffold_chain.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * LIS element (Longest Increasing Subsequence)
 */
typedef struct {
  int coverage;
  int sparseness;
  int next_region;
} match_scaffold_lis_t;
int match_scaffold_lis_cmp(
    const match_scaffold_lis_t* const a,
    const match_scaffold_lis_t* const b) {
  int coverage_diff = b->coverage - a->coverage;
  if (coverage_diff) return coverage_diff;
  return a->sparseness - b->sparseness;
}
/*
 * Compatible region
 */
bool match_scaffold_chain_is_compatible(
    match_alignment_region_t* const region_a,
    match_alignment_region_t* const region_b) {
  if (match_alignment_region_text_overlap(region_a,region_b)) return false;
  //if (match_alignment_region_key_overlap(current_region,next_region) return false; // Guaranteed
  return (match_alignment_region_cmp_key_position(region_a,region_b) < 0);
}
/*
 * Compute LIS
 */
uint64_t match_scaffold_chain_compute_lis(
    match_scaffold_t* const match_scaffold,
    match_scaffold_lis_t* const lis_vector,
    const bool left_gap_alignment) {
  // Parameters
  match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions;
  const int64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  // Compute the LIS-vector
  int64_t region_idx, best_lis_idx = num_alignment_regions-1;
  for (region_idx=num_alignment_regions-1;region_idx>=0;--region_idx) {
    match_alignment_region_t* const current_region = match_alignment_region + region_idx;
    const uint64_t current_region_coverage = match_alignment_region_text_coverage(current_region);
    match_scaffold_lis_t best_partial_lis = { // LIS containing only current region
        .coverage = current_region_coverage,
        .sparseness = 0,
        .next_region = num_alignment_regions, // No-Next
    };
    // Search for the best compatible chain-region (non-overlapping & in-order)
    int64_t closest_reachable_idx = num_alignment_regions;
    int64_t next_region_idx = region_idx + 1;
    while (next_region_idx < closest_reachable_idx) {
      // Test compatible
      match_alignment_region_t* const next_region = match_alignment_region + next_region_idx;
      if (match_scaffold_chain_is_compatible(current_region,next_region)) {
        // Compare with chained LIS
        const match_scaffold_lis_t chained_lis = {
            .coverage = current_region_coverage + lis_vector[next_region_idx].coverage,
            .sparseness = match_alignment_region_text_distance(current_region,next_region) +
                          lis_vector[next_region_idx].sparseness,
            .next_region = next_region_idx,
        };
        int lis_cmp = match_scaffold_lis_cmp(&chained_lis,&best_partial_lis) < 0;
        if (lis_cmp || (left_gap_alignment && lis_cmp==0)) {
          best_partial_lis = chained_lis;
          closest_reachable_idx = MIN(closest_reachable_idx,lis_vector[next_region_idx].next_region);
        }
      }
      ++next_region_idx;
    }
    // Set current LIS
    lis_vector[region_idx] = best_partial_lis;
    // Compare with the best LIS so far
    int lis_cmp = match_scaffold_lis_cmp(lis_vector+region_idx,lis_vector+best_lis_idx);
    if (lis_cmp || (left_gap_alignment && lis_cmp==0)) {
      best_lis_idx = region_idx;
    }
  }
  // Return index of the best LIS
  return best_lis_idx;
}
void match_scaffold_chain_backtrace_lis(
    match_scaffold_t* const match_scaffold,
    match_scaffold_lis_t* const lis_vector,
    const uint64_t best_lis_idx) {
  // Parameters
  match_alignment_region_t* match_alignment_region = match_scaffold->alignment_regions;
  const int64_t lis_vector_length = match_scaffold->num_alignment_regions;
  // Traverse the LIS and store the best chain of regions
  uint64_t region_idx = 0, i = best_lis_idx;
  while (i < lis_vector_length) {
    match_alignment_region[region_idx] = match_alignment_region[i];
    ++region_idx;
    i = lis_vector[i].next_region;
  }
  match_scaffold->num_alignment_regions = region_idx;
  match_scaffold->scaffolding_coverage = lis_vector[best_lis_idx].coverage;
}
void match_scaffold_chain_alignment_regions(
    match_scaffold_t* const match_scaffold,
    const bool left_gap_alignment,
    mm_allocator_t* const mm_allocator) {
  PROFILE_START(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
  // Sort alignment regions by text-offsets
  match_scaffold_sort_alignment_regions(match_scaffold);
  // Allocate DP-LIS table
  mm_allocator_push_state(mm_allocator);
  const int64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  match_scaffold_lis_t* const lis_vector = mm_allocator_calloc(
      mm_allocator,num_alignment_regions,match_scaffold_lis_t,true);
  // Compute LIS (longest increasing sequence of alignment-regions)
  const uint64_t best_lis_idx =
      match_scaffold_chain_compute_lis(
          match_scaffold,lis_vector,left_gap_alignment);
  // Keep the best chain (given by the LIS)
  match_scaffold_chain_backtrace_lis(match_scaffold,lis_vector,best_lis_idx);
  mm_allocator_pop_state(mm_allocator);
}
/*
 * Exact extend alignment-regions
 */
uint64_t match_scaffold_chain_exact_extend_right(
    match_scaffold_t* const match_scaffold,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t region_idx) {
  match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + region_idx;
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  uint64_t inc_coverage = 0;
  int64_t right_key_max, right_text_max;
  if (region_idx==num_alignment_regions-1) {
    right_key_max = key_length-1;
    right_text_max = text_length-1;
  } else {
    match_alignment_region_t* const next_alignment_region = match_scaffold->alignment_regions + (region_idx+1);
    right_key_max = match_alignment_region_get_key_begin(next_alignment_region)-1;
    right_text_max = match_alignment_region_get_text_begin(next_alignment_region)-1;
  }
  int64_t right_key = match_alignment_region_get_key_end(match_alignment_region);
  int64_t right_text = match_alignment_region_get_text_end(match_alignment_region);
  while (right_key_max>=right_key && right_text_max>=right_text) {
    // Check match
    const uint8_t candidate_enc = text[right_text];
    if (candidate_enc == ENC_DNA_CHAR_N || candidate_enc != key[right_key]) break;
    ++right_key;
    ++right_text;
    ++inc_coverage;
  }
  match_alignment_region_set_key_end(match_alignment_region,right_key);
  match_alignment_region_set_text_end(match_alignment_region,right_text);
  return inc_coverage;
}
uint64_t match_scaffold_chain_exact_extend_left(
    match_scaffold_t* const match_scaffold,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t region_idx) {
  match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + region_idx;
  const uint64_t region_key_begin = match_alignment_region_get_key_begin(match_alignment_region);
  const uint64_t region_text_begin = match_alignment_region_get_text_begin(match_alignment_region);
  uint64_t inc_coverage = 0;
  if (region_key_begin>0 && region_text_begin>0) {
    int64_t left_key_max, left_text_max;
    if (region_idx==0) {
      left_key_max = 0;
      left_text_max = 0;
    } else {
      match_alignment_region_t* const prev_alignment_region = match_scaffold->alignment_regions + (region_idx-1);
      left_key_max = match_alignment_region_get_key_end(prev_alignment_region);
      left_text_max = match_alignment_region_get_text_end(prev_alignment_region);
    }
    int64_t left_key = region_key_begin-1;
    int64_t left_text = region_text_begin-1;
    while (left_key_max<=left_key && left_text_max<=left_text) {
      // Check match
      const uint8_t candidate_enc = text[left_text];
      if (candidate_enc == ENC_DNA_CHAR_N || candidate_enc != key[left_key]) break;
      --left_key;
      --left_text;
      ++inc_coverage;
    }
    match_alignment_region_set_key_begin(match_alignment_region,left_key+1);
    match_alignment_region_set_text_begin(match_alignment_region,left_text+1);
  }
  return inc_coverage;
}
void match_scaffold_chain_exact_extend(
    match_scaffold_t* const match_scaffold,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const bool left_gap_alignment) {
  // Parameters
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  uint64_t inc_coverage = 0;
  int64_t i;
  // Switch gap alignment
  if (left_gap_alignment) {
    // Exact extend all alignment-regions (from right to left)
    for (i=num_alignment_regions-1;i>=0;--i) {
      // Extend right
      inc_coverage += match_scaffold_chain_exact_extend_right(
          match_scaffold,key,key_length,text,text_length,i);
      // Extend left
      inc_coverage += match_scaffold_chain_exact_extend_left(
          match_scaffold,key,key_length,text,text_length,i);
    }
  } else {
    // Exact extend all alignment-regions (from left to right)
    for (i=0;i<num_alignment_regions;++i) {
      // Extend left
      inc_coverage += match_scaffold_chain_exact_extend_left(
          match_scaffold,key,key_length,text,text_length,i);
      // Extend right
      inc_coverage += match_scaffold_chain_exact_extend_right(
          match_scaffold,key,key_length,text,text_length,i);
    }
  }
  match_scaffold->scaffolding_coverage += inc_coverage;
}
/*
 * Alignment-region chain scaffolding
 */
void match_scaffold_chain(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool exact_extend,
    mm_allocator_t* const mm_allocator) {
  PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_SCAFFOLDS);
  PROFILE_START(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
  // Parameters
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->text_length;
  // Find a compatible chain of alignment-regions
  if (match_scaffold->num_alignment_regions > 0) {
    const bool left_gap_alignment = (text_trace->strand==Forward);
    match_scaffold_chain_alignment_regions(match_scaffold,left_gap_alignment,mm_allocator);
    if (match_scaffold->num_alignment_regions > 0) {
      // Extend alignment-regions as to maximize coverage
      if (exact_extend && match_scaffold->scaffolding_coverage < key_length) {
        match_scaffold_chain_exact_extend(
            match_scaffold,key,key_length,
            text,text_length,left_gap_alignment);
      }
      match_scaffold->match_alignment.score =
          key_length - match_scaffold->scaffolding_coverage; // Set score as matching bases
      PROFILE_STOP(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
      return;
    }
  }
  // Set score & coverage
  match_scaffold->match_alignment.score = key_length;
  match_scaffold->scaffolding_coverage = 0;
  PROFILE_STOP(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
}
