/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_region_chain.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_scaffold_region_chain.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Exact extend matching regions
 */
void match_scaffold_exact_extend(
    matches_t* const matches,match_scaffold_t* const match_scaffold,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const uint64_t text_length,
    const bool* const allowed_enc) {
  PROFILE_START(GP_MATCH_SCAFFOLD_EXTEND_REGIONS,PROFILE_LEVEL);
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
  PROFILE_STOP(GP_MATCH_SCAFFOLD_EXTEND_REGIONS,PROFILE_LEVEL);
}
/*
 * Chain matching regions
 *   @match_scaffold->num_scaffold_regions
 *   @match_scaffold->scaffold_regions
 */
void match_scaffold_chain_matching_regions(
    matches_t* const matches,match_scaffold_t* const match_scaffold,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const bool* const allowed_enc,
    const uint64_t max_error,mm_stack_t* const stack) {
  PROFILE_START(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  // Sort matching regions
  match_scaffold_sort_regions_matching(match_scaffold);
  // Check overlapping
  bool overlapping_text = false, unsorted_read = false, max_diff = false;
  region_matching_t* last_region_matching = NULL;
  uint64_t i, coverage = 0;
  for (i=0;i<num_scaffold_regions;++i) {
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    if (last_region_matching!=NULL) {
      if (last_region_matching->text_end > region_matching->text_begin) {
        overlapping_text = true; break; // Might be an deletion in the reference
      }
      if (last_region_matching->key_end > region_matching->key_begin) {
        unsorted_read = true; break;  // Might be a deletion in the read
      }
      const int64_t text_gap = region_matching->text_begin - last_region_matching->text_end;
      const int64_t key_gap = region_matching->key_begin - last_region_matching->key_end;
      const int64_t diff = text_gap-key_gap;
      if (ABS(diff) > max_error) {
        max_diff = true; break;
      }
    }
    coverage += region_matching->key_end - region_matching->key_begin;
    last_region_matching = region_matching;
  }
  if (overlapping_text || unsorted_read || max_diff) {
    match_scaffold->num_scaffold_regions = 0; // Disable region chaining
    match_scaffold->scaffolding_coverage = 0;
  } else {
    match_scaffold->scaffolding_coverage = coverage;
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
}
/*
 * Region-Chain Scaffolding
 */
bool match_scaffold_region_chain(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  // Init
  match_scaffold->scaffold_type = scaffold_region_chain;
  // Check number of regions to chain/extend
  if (match_scaffold->num_scaffold_regions > 0) {
    // Parameters
    const uint8_t* const key = align_input->key;
    const uint64_t key_length = align_input->key_length;
    const uint8_t* const text = align_input->text;
    const uint64_t text_length = align_input->text_length;
    const uint64_t max_error = align_parameters->max_error;
    const bool* const allowed_enc = align_parameters->allowed_enc;
    // Find a compatible chain of matching-regions
    match_scaffold_chain_matching_regions(matches,match_scaffold,key,key_length,text,allowed_enc,max_error,mm_stack);
    // Extend matching-regions as to maximize coverage
    if (match_scaffold->num_scaffold_regions > 0) {
      PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_SUCCESS);
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE,(100*match_scaffold->scaffolding_coverage)/key_length);
      if (match_scaffold->scaffolding_coverage < key_length) {
        match_scaffold_exact_extend(matches,match_scaffold,key,key_length,text,text_length,allowed_enc);
      }
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EXTEND_REGIONS_COVERAGE,(100*match_scaffold->scaffolding_coverage)/key_length);
      // Set score as matching bases
      match_scaffold->match_alignment.score = key_length - match_scaffold->scaffolding_coverage;
      // Return OK
      return true;
    }
  }
  // Set score & coverage
  match_scaffold->match_alignment.score = align_input->key_length;
  match_scaffold->scaffolding_coverage = 0;
  // Return fail
  return false;
}
