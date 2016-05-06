/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/match_scaffold.h"
#include "matches/match_scaffold_region_chain.h"
#include "matches/match_scaffold_levenshtein.h"
#include "matches/matches_cigar.h"
#include "matches/match_align_rl.h"
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
  match_scaffold->scaffold_regions_rl = false;
  // Scaffold Matching Region
  match_scaffold->scaffold_regions = NULL;
  match_scaffold->num_scaffold_regions = 0;
  match_scaffold->scaffolding_coverage = 0;
  // Underlying Alignment
  match_scaffold->match_alignment.cigar_length = 0;
}
/*
 * Accessors
 */
bool match_scaffold_is_null(match_scaffold_t* const match_scaffold) {
  return match_scaffold->num_scaffold_regions==0;
}
/*
 * RL-Translation
 */
void match_scaffold_rl_translate_regions(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    matches_t* const matches) {
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  uint64_t i;
  for (i=0;i<num_scaffold_regions;++i) {
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    // Translate into Text-Space
    region_matching->matching_type = region_matching_approximate;
    // Translate CIGAR
    match_align_rl_translate_region_cigar(region_matching,align_input,
        align_parameters->left_gap_alignment,matches->cigar_vector);
    // Translate offsets
    uint32_t* const rl_key_runs_acc = align_input->rl_key_runs_acc;
    region_matching->key_begin = archive_text_rl_get_decoded_offset_exl(rl_key_runs_acc,region_matching->key_begin);
    region_matching->key_end = archive_text_rl_get_decoded_offset_exl(rl_key_runs_acc,region_matching->key_end);
    uint32_t* const rl_text_runs_acc = align_input->rl_text_runs_acc;
    region_matching->text_begin = archive_text_rl_get_decoded_offset_exl(rl_text_runs_acc,region_matching->text_begin);
    region_matching->text_end = archive_text_rl_get_decoded_offset_exl(rl_text_runs_acc,region_matching->text_end);
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
 *   @match_scaffold->num_scaffold_regions
 *   @match_scaffold->scaffold_regions
 */
void match_scaffold_adaptive(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    matches_t* const matches,
    mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
  // Translate matching regions
  const bool rl_space = match_scaffold->scaffold_regions_rl;
  if (rl_space) {
    match_scaffold_rl_translate_regions(match_scaffold,align_input,align_parameters,matches);
    match_scaffold->scaffold_regions_rl = false;
  }
  // Select proper scaffold approach
  switch (match_scaffold->scaffold_type) {
    case scaffold_none:
      // Scaffold chaining matching regions (from region-profile)
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_MATCHING_REGIONS,match_scaffold->num_scaffold_regions);
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_MATCHING_COVERAGE,
          (100*match_scaffold->scaffolding_coverage)/align_input->key_length);
      match_scaffold->scaffold_type = scaffold_region_chain;
      match_scaffold_region_chain(match_scaffold,align_input,align_parameters,!rl_space,mm_stack);
      // Check coverage
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE,
          (100*match_scaffold->scaffolding_coverage)/align_input->key_length);
      const uint64_t max_coverage_bound = BOUNDED_SUBTRACTION(
          align_input->key_length,align_input->region_alignment->distance_min_bound,0);
      if (max_coverage_bound >= align_parameters->global_min_identity &&
          match_scaffold->scaffolding_coverage >= max_coverage_bound) break;
      // no break
    case scaffold_region_chain:
      // Scaffold from Levenshtein-alignment
      match_scaffold->scaffold_type = scaffold_levenshtein;
      match_scaffold_levenshtein(match_scaffold,align_input,align_parameters,matches,mm_stack);
      match_scaffold_region_chain(match_scaffold,align_input,align_parameters,false,mm_stack);
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
int region_matching_cmp_text_position(const region_matching_t* const a,const region_matching_t* const b) {
  return a->text_begin - b->text_begin;
}
void match_scaffold_sort_regions_matching(match_scaffold_t* const match_scaffold) {
  // Sort Scaffold regions (region_matching_t) wrt their starting position in the text
  void* array = match_scaffold->scaffold_regions;
  const size_t count = match_scaffold->num_scaffold_regions;
  qsort(array,count,sizeof(region_matching_t),(int (*)(const void *,const void *))region_matching_cmp_text_position);
}
/*
 * Display
 */
void match_scaffold_print_matching_region(
    FILE* const stream,matches_t* const matches,
    const uint64_t region_matching_id,
    region_matching_t* const region_matching) {
  // Print matching region
  switch (region_matching->matching_type) {
    case region_matching_exact: tab_fprintf(stream,"    %"PRIu64"[exact]\t",region_matching_id); break;
    case region_matching_approximate: tab_fprintf(stream,"    %"PRIu64"[approximate]\t",region_matching_id); break;
    default: GEM_INVALID_CASE(); break;
  }
  tab_fprintf(stream,"-> [%"PRIu64",%"PRIu64") ~> [+%"PRIu64",+%"PRIu64")",
      region_matching->key_begin,region_matching->key_end,
      region_matching->text_begin,region_matching->text_end);
  // Print CIGAR
  if (matches!=NULL && region_matching->matching_type==region_matching_approximate && region_matching->cigar_length>0) {
    tab_fprintf(stream,"\tCIGAR=");
    match_cigar_print(stream,matches->cigar_vector,region_matching->cigar_buffer_offset,region_matching->cigar_length);
  }
  tab_fprintf(stream,"\n");
}
void match_scaffold_print(
    FILE* const stream,
    matches_t* const matches,
    match_scaffold_t* const match_scaffold) {
  tab_fprintf(stream,"[GEM]>Matching.Scaffold.Regions\n");
  switch (match_scaffold->scaffold_type) {
    case scaffold_none: tab_fprintf(stream,"  => Scaffold.type -None-\n"); break;
    case scaffold_region_chain: tab_fprintf(stream,"  => Scaffold.type region.chain\n"); break;
    case scaffold_levenshtein: tab_fprintf(stream,"  => Scaffold.type levenshtein\n"); break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  tab_fprintf(stream,"  => Num.scaffold.regions %"PRIu64"\n",num_scaffold_regions);
  tab_fprintf(stream,"  => Scaffold.coverage %"PRIu64"\n",match_scaffold->scaffolding_coverage);
  uint64_t i;
  for (i=0;i<num_scaffold_regions;++i) {
    // Print matching region
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    match_scaffold_print_matching_region(stream,matches,i,region_matching);
  }
}
void match_scaffold_print_pretty_matching_region(
    FILE* const stream,
    region_matching_t* const region_matching,
    uint8_t* key,
    uint64_t key_length,
    uint8_t* text) {
  // Compute offset
  uint64_t i, offset_text = 0, offset_key = 0;
  if (region_matching->key_begin > region_matching->text_begin) {
    offset_key = region_matching->key_begin-region_matching->text_begin;
  } else {
    offset_text = region_matching->text_begin-region_matching->key_begin;
    for (i=0;i<offset_text;++i) fprintf(stream," ");
  }
  // Print Key
  for (i=offset_key;i<key_length;++i) {
    if (region_matching->key_begin <= i && i < region_matching->key_end) {
      if (text[offset_text+i]==key[i]) {
        fprintf(stream,"%c",dna_decode(key[i]));
      } else {
        fprintf(stream,"%c",tolower(dna_decode(key[i])));
      }
    } else {
      fprintf(stream,"-");
    }
  }
  fprintf(stream,"\n");
}
void match_scaffold_print_pretty(
    FILE* const stream,
    matches_t* const matches,
    match_scaffold_t* const match_scaffold,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    mm_stack_t* const mm_stack) {
  tab_fprintf(stream,"[GEM]>Match.Scaffold\n");
  // Print regions matching
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  uint64_t i;
  for (i=0;i<num_scaffold_regions;++i) {
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    match_scaffold_print_matching_region(stream,matches,i,region_matching);
    match_scaffold_print_pretty_matching_region(stream,region_matching,key,key_length,text);
  }
  // Print text
  dna_buffer_print(stream,text,text_length,false);
  fprintf(stream,"\n");
}

