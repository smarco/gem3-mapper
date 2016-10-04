/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/align/match_alignment_region_rl.h"
#include "matches/scaffold/match_scaffold_chain.h"
#include "matches/scaffold/match_scaffold.h"
#include "matches/scaffold/match_scaffold_levenshtein.h"
#include "matches/matches_cigar.h"
#include "archive/archive_text_rl.h"

/*
 * Debug
 */
#define DEBUG_MATCH_SCAFFOLD       GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Setup
 */
void match_scaffold_init(match_scaffold_t* const match_scaffold) {
  // Scaffold Properties
  match_scaffold->scaffold_type = scaffold_none;
  match_scaffold->alignment_regions_rl = false;
  // Scaffold Alignment-Regions
  match_scaffold->alignment_regions = NULL;
  match_scaffold->num_alignment_regions = 0;
  match_scaffold->scaffolding_coverage = 0;
  // Underlying Alignment
  match_scaffold->match_alignment.cigar_length = 0;
}
/*
 * Accessors
 */
bool match_scaffold_is_null(match_scaffold_t* const match_scaffold) {
  return match_scaffold->num_alignment_regions==0;
}
/*
 * RL-Translation
 */
void match_scaffold_rl_translate_regions(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    matches_t* const matches) {
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  uint64_t i;
  for (i=0;i<num_alignment_regions;++i) {
    match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + i;
    // Translate into Text-Space
    match_alignment_region->region_type = match_alignment_region_approximate;
    // Translate CIGAR
    match_alignment_region_rl_translate(match_alignment_region,align_input,
        align_parameters->left_gap_alignment,matches->cigar_vector);
    // Translate offsets
    uint32_t* const rl_key_runs_acc = align_input->rl_key_runs_acc;
    match_alignment_region->key_begin = archive_text_rl_get_decoded_offset_exl(rl_key_runs_acc,match_alignment_region->key_begin);
    match_alignment_region->key_end = archive_text_rl_get_decoded_offset_exl(rl_key_runs_acc,match_alignment_region->key_end);
    uint32_t* const rl_text_runs_acc = align_input->rl_text_runs_acc;
    match_alignment_region->text_begin = archive_text_rl_get_decoded_offset_exl(rl_text_runs_acc,match_alignment_region->text_begin);
    match_alignment_region->text_end = archive_text_rl_get_decoded_offset_exl(rl_text_runs_acc,match_alignment_region->text_end);
  }
}
/*
 * Adaptive Scaffolding of the alignment (Best effort)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->bpm_pattern_tiles
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->key_trim_left
 *   @align_input->key_trim_right
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_input->text_offset_base_begin
 *   @align_input->text_offset_base_end
 *   @align_parameters->max_error
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->scaffolding_min_coverage
 *   @align_parameters->scaffolding_matching_min_length
 *   @align_parameters->scaffolding_homopolymer_min_context
 *   @match_scaffold->num_alignment_regions
 *   @match_scaffold->alignment_regions
 */
void match_scaffold_adaptive(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    matches_t* const matches,
    mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
  // Translate alignment-regions
  const bool rl_space = match_scaffold->alignment_regions_rl;
  if (rl_space) {
    match_scaffold_rl_translate_regions(match_scaffold,align_input,align_parameters,matches);
    match_scaffold->alignment_regions_rl = false;
  }
  // Select proper scaffold approach
  switch (match_scaffold->scaffold_type) {
    case scaffold_none:
      // Scaffold chaining alignment-regions (from region-profile)
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_REGIONS,match_scaffold->num_alignment_regions);
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_COVERAGE,
          (100*match_scaffold->scaffolding_coverage)/align_input->key_length);
      match_scaffold->scaffold_type = scaffold_alignment_chain;
      match_scaffold_chain(match_scaffold,align_input,align_parameters,!rl_space,mm_stack);
      // Check coverage
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE,
          (100*match_scaffold->scaffolding_coverage)/align_input->key_length);
      const uint64_t max_coverage_bound = BOUNDED_SUBTRACTION(
          align_input->key_length,align_input->alignment->distance_min_bound,0);
      if (max_coverage_bound >= align_parameters->global_min_identity &&
          match_scaffold->scaffolding_coverage >= max_coverage_bound) break;
      // no break
    case scaffold_alignment_chain:
      // Scaffold from Levenshtein-alignment
      match_scaffold->scaffold_type = scaffold_levenshtein;
      match_scaffold_levenshtein(match_scaffold,align_input,align_parameters,matches,mm_stack);
      match_scaffold_chain(match_scaffold,align_input,align_parameters,false,mm_stack);
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_COVERAGE,
          (100*match_scaffold->scaffolding_coverage)/align_input->key_length);
      // no break
    default:
      break;
  }
  //  // DEBUG
  //  match_scaffold_print_pretty(stderr,matches,match_scaffold,
  //      align_input->key,align_input->key_length,
  //      align_input->text,align_input->text_length,mm_stack);
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold (scaffold alignment)\n");
    tab_global_inc();
    match_scaffold_print(stderr,matches,match_scaffold);
    tab_global_dec();
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
}
/*
 * Sorting
 */
#define VECTOR_SORT_NAME                 match_alignment_region
#define VECTOR_SORT_TYPE                 match_alignment_region_t
#define VECTOR_SORT_CMP(a,b)             match_alignment_region_cmp_text_position(a,b)
#include "utils/vector_sort.h"
void match_scaffold_sort_alignment_regions(match_scaffold_t* const match_scaffold) {
  buffer_sort_match_alignment_region(match_scaffold->alignment_regions,match_scaffold->num_alignment_regions);
}
/*
 * Display
 */
void match_scaffold_print(
    FILE* const stream,
    matches_t* const matches,
    match_scaffold_t* const match_scaffold) {
  tab_fprintf(stream,"[GEM]>Scaffold.Alignment.Regions\n");
  switch (match_scaffold->scaffold_type) {
    case scaffold_none: tab_fprintf(stream,"  => Scaffold.type -None-\n"); break;
    case scaffold_alignment_chain: tab_fprintf(stream,"  => Scaffold.type alignment.chain\n"); break;
    case scaffold_levenshtein: tab_fprintf(stream,"  => Scaffold.type levenshtein\n"); break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  tab_fprintf(stream,"  => Num.scaffold.alignment.regions %"PRIu64"\n",num_alignment_regions);
  tab_fprintf(stream,"  => Scaffold.coverage              %"PRIu64"\n",match_scaffold->scaffolding_coverage);
  uint64_t i;
  for (i=0;i<num_alignment_regions;++i) {
    // Print alignment-regions
    match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + i;
    match_alignment_region_print(stream,match_alignment_region,i,matches);
  }
}
void match_scaffold_print_pretty(
    FILE* const stream,
    matches_t* const matches,
    match_scaffold_t* const match_scaffold,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length) {
  tab_fprintf(stream,"[GEM]>Match.Scaffold\n");
  // Print alignment-regions
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  uint64_t i;
  for (i=0;i<num_alignment_regions;++i) {
    match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + i;
    match_alignment_region_print(stream,match_alignment_region,i,matches);
    match_alignment_region_print_pretty(stream,match_alignment_region,key,key_length,text);
  }
  // Print text
  dna_buffer_print(stream,text,text_length,false);
  fprintf(stream,"\n");
}

