/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_scaffold.h"
#include "match_scaffold_region_chain.h"
#include "match_scaffold_levenshtein.h"
#include "matches_cigar.h"

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
 * Adaptive Scaffolding of the alignment (Best effort)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->max_error
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->scaffolding_min_coverage
 *   @align_parameters->scaffolding_matching_min_length
 *   @align_parameters->scaffolding_homopolymer_min_context
 *   @match_scaffold->num_scaffold_regions
 *   @match_scaffold->scaffold_regions
 */
void match_scaffold_adaptive(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
  // Select proper scaffold approach
  switch (match_scaffold->scaffold_type) {
    case scaffold_none:
      // Chaining Matching Regions
      match_scaffold_region_chain(matches,align_input,align_parameters,match_scaffold,mm_stack);
      // no break
    case scaffold_region_chain:
      // Scaffold from Levenshtein-alignment
      if (match_scaffold->scaffolding_coverage < align_parameters->scaffolding_min_coverage) {
        match_scaffold_levenshtein(matches,align_input,align_parameters,match_scaffold,mm_stack);
      }
      // no break
    default: break;
  }
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
void match_scaffold_print(
    FILE* const stream,matches_t* const matches,match_scaffold_t* const match_scaffold) {
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  tab_fprintf(stream,"[GEM]>Matching.Scaffolded.Regions\n");
  switch (match_scaffold->scaffold_type) {
    case scaffold_none: tab_fprintf(stream,"  => Scaffold.type -None-\n"); return;
    case scaffold_region_chain: tab_fprintf(stream,"  => Scaffold.type region.chain\n"); break;
    case scaffold_levenshtein: tab_fprintf(stream,"  => Scaffold.type levenshtein\n"); break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  tab_fprintf(stream,"  => Num.scaffold.regions %"PRIu64"\n",num_scaffold_regions);
  tab_fprintf(stream,"  => Scaffold.coverage %"PRIu64"\n",match_scaffold->scaffolding_coverage);
  uint64_t i;
  for (i=0;i<num_scaffold_regions;++i) {
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    // Print matching region
    switch (region_matching->matching_type) {
      case region_matching_exact: tab_fprintf(stream,"    %"PRIu64"[exact]\t",i); break;
      case region_matching_approximate: tab_fprintf(stream,"    %"PRIu64"[approximate]\t",i); break;
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
}

