/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_region_chain.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/match_scaffold_region_chain.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * LIS element (Longest Increasing Subsequence)
 */
typedef struct {
  uint64_t coverage;
  uint64_t sparseness;
  uint64_t next_region;
  bool region_chained;
} match_scaffold_lis_t;
int match_scaffold_lis_cmp(
    const match_scaffold_lis_t* const a,
    const match_scaffold_lis_t* const b) {
  int coverage_diff =  b->coverage - a->coverage;
  if (coverage_diff) return coverage_diff;
  return a->sparseness - b->sparseness;
}
/*
 * Compute LIS
 */
void match_scaffold_chain_compute_lis(
    match_scaffold_t* const match_scaffold,
    match_scaffold_lis_t* const lis_vector) {
  // Parameters
  region_matching_t* const region_matching = match_scaffold->scaffold_regions;
  const int64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  // Compute the LIS-vector
  int64_t region_idx;
  for (region_idx=num_scaffold_regions-1;region_idx>=0;--region_idx) {
    // Init
    region_matching_t* current_region = region_matching + region_idx;
    const uint64_t current_region_coverage = region_matching_text_coverage(current_region);
    match_scaffold_lis_t best_lis = { // LIS containing only current region
        .coverage = current_region_coverage,
        .sparseness = 0,
        .next_region = num_scaffold_regions,
        .region_chained = true,
    };
    // Search for the first compatible region (non-overlapping & in-order)
    int64_t next_region_idx = region_idx + 1;
    region_matching_t* next_region = current_region + 1;
    while (next_region_idx < num_scaffold_regions) {
      // Skip unchained regions
      if (lis_vector[next_region_idx].region_chained) {
        // Test region compatible
        if (!region_matching_text_overlap(current_region,next_region) &&
            region_matching_key_cmp(current_region,next_region) < 0) {
          break;
        }
        // Pick best LIS (from discarding current region)
        if (match_scaffold_lis_cmp(&lis_vector[next_region_idx],&best_lis) < 0) {
          best_lis.coverage = lis_vector[next_region_idx].coverage;
          best_lis.sparseness = lis_vector[next_region_idx].sparseness;
          best_lis.next_region = next_region_idx;
          best_lis.region_chained = false; // Discards current region
        }
      }
      ++next_region_idx;
      ++next_region;
    }
    // Next region compatible with current
    if (next_region_idx < num_scaffold_regions) {
      // Join current-region with first non-overlapping
      const match_scaffold_lis_t joint_lis = {
          .coverage = current_region_coverage + lis_vector[next_region_idx].coverage,
          .sparseness = region_matching_text_distance(current_region,next_region) +
                        lis_vector[next_region_idx].sparseness,
          .next_region = next_region_idx,
          .region_chained = true,
      };
      // Pick best LIS
      if (match_scaffold_lis_cmp(&joint_lis,&best_lis) < 0) {
        best_lis = joint_lis;
      }
    }
    // Set current LIS
    lis_vector[region_idx] = best_lis;
  }
}
void match_scaffold_chain_store_lis(
    match_scaffold_t* const match_scaffold,
    match_scaffold_lis_t* const lis_vector) {
  // Parameters
  region_matching_t* region_matching = match_scaffold->scaffold_regions;
  const int64_t lis_vector_length = match_scaffold->num_scaffold_regions;
  // Traverse the LIS and store the best chain of regions
  uint64_t num_scaffold_regions = 0, i = 0;
  while (i < lis_vector_length) {
    if (lis_vector[i].region_chained) {
      region_matching[num_scaffold_regions] = region_matching[i];
      ++num_scaffold_regions;
    }
    i = lis_vector[i].next_region;
  }
  match_scaffold->num_scaffold_regions = num_scaffold_regions;
  match_scaffold->scaffolding_coverage = lis_vector[0].coverage;
}
void match_scaffold_chain_matching_regions(
    match_scaffold_t* const match_scaffold,
    mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
  // Sort matching regions by text-offsets
  match_scaffold_sort_regions_matching(match_scaffold);
  // Allocate DP-LIS table
  mm_stack_push_state(mm_stack);
  const int64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  match_scaffold_lis_t* const lis_vector = mm_stack_calloc(
      mm_stack,num_scaffold_regions,match_scaffold_lis_t,true);
  // Compute LIS (longest increasing sequence of matching-regions)
  match_scaffold_chain_compute_lis(match_scaffold,lis_vector);
  // Keep the best chain (given by the LIS)
  match_scaffold_chain_store_lis(match_scaffold,lis_vector);
  mm_stack_pop_state(mm_stack);

}
/*
 * Exact extend matching regions
 */
void match_scaffold_exact_extend(
    match_scaffold_t* const match_scaffold,
    const bool* const allowed_enc,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length) {
  // Extend all matching regions (Exact extend of matching regions)
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  const uint64_t last_region = num_scaffold_regions-1;
  uint64_t i, inc_coverage = 0;
  for (i=0;i<num_scaffold_regions;++i) {
    // Try to left extend
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    const int64_t left_key_max = (i==0) ? 0 : match_scaffold->scaffold_regions[i-1].key_end;
    const int64_t left_text_max = (i==0) ? 0 : match_scaffold->scaffold_regions[i-1].text_end;
    int64_t left_key = region_matching->key_begin-1;
    int64_t left_text = region_matching->text_begin-1;
    while (left_key_max<=left_key && left_text_max<=left_text) {
      // Check match
      const uint8_t candidate_enc = text[left_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[left_key]) break;
      --left_key;
      --left_text;
      ++inc_coverage;
    }
    region_matching->key_begin = left_key+1;
    region_matching->text_begin = left_text+1;
    // Try to right extend
    const int64_t right_key_max = (i==last_region) ? key_length-1 : match_scaffold->scaffold_regions[i+1].key_begin-1;
    const int64_t right_text_max = (i==last_region) ? text_length-1 : match_scaffold->scaffold_regions[i+1].text_begin-1;
    int64_t right_key = region_matching->key_end;
    int64_t right_text = region_matching->text_end;
    while (right_key_max>=right_key && right_text_max>=right_text) {
      // Check match
      const uint8_t candidate_enc = text[right_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[right_key]) break;
      ++right_key;
      ++right_text;
      ++inc_coverage;
    }
    region_matching->key_end = right_key;
    region_matching->text_end = right_text;
  }
  match_scaffold->scaffolding_coverage += inc_coverage;
}

/*
 * Region-Chain Scaffolding
 */
void match_scaffold_region_chain(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool exact_extend,
    mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_SCAFFOLDS);
  PROFILE_START(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  const uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  const bool* const allowed_enc = align_parameters->allowed_enc;
  // Find a compatible chain of matching-regions
  if (match_scaffold->num_scaffold_regions > 0) {
    match_scaffold_chain_matching_regions(match_scaffold,mm_stack);
    if (match_scaffold->num_scaffold_regions > 0) {
      // Extend matching-regions as to maximize coverage
      if (exact_extend && match_scaffold->scaffolding_coverage < key_length) {
        match_scaffold_exact_extend(match_scaffold,allowed_enc,key,key_length,text,text_length);
      }
      match_scaffold->match_alignment.score =
          key_length - match_scaffold->scaffolding_coverage; // Set score as matching bases
      PROFILE_STOP(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
      return;
    }
  }
  // Set score & coverage
  match_scaffold->match_alignment.score = align_input->key_length;
  match_scaffold->scaffolding_coverage = 0;
  PROFILE_STOP(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
}
