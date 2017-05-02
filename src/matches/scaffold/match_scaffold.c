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

#include "matches/align/match_alignment_region_rl.h"
#include "matches/scaffold/match_scaffold_chain.h"
#include "matches/scaffold/match_scaffold.h"
#include "matches/scaffold/match_scaffold_levenshtein.h"
#include "matches/matches_cigar.h"
#include "archive/archive_text_rl.h"

/*
 * Debug
 */
#define DEBUG_MATCH_SCAFFOLD        GEM_DEEP_DEBUG
#define DEBUG_MATCH_SCAFFOLD_CHECK  false

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Setup
 */
void match_scaffold_init(
    match_scaffold_t* const match_scaffold) {
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
void match_scaffold_destroy(
    match_scaffold_t* const match_scaffold,
    mm_allocator_t* const mm_allocator) {
  if (match_scaffold->alignment_regions!=NULL) {
    match_scaffold_free_alignment_region(
        match_scaffold,match_scaffold->alignment_regions,mm_allocator);
  }
}
/*
 * Accessors
 */
bool match_scaffold_is_null(match_scaffold_t* const match_scaffold) {
  return match_scaffold->num_alignment_regions==0;
}
match_alignment_region_t* match_scaffold_allocate_alignment_region(
    match_scaffold_t* const match_scaffold,
    const uint64_t num_alignment_regions,
    mm_allocator_t* const mm_allocator) {
  return mm_allocator_malloc(mm_allocator,num_alignment_regions*sizeof(match_alignment_region_t));
}
void match_scaffold_free_alignment_region(
    match_scaffold_t* const match_scaffold,
    match_alignment_region_t* const match_alignment_region,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,match_alignment_region);
}
/*
 * RL-Translation
 */
void match_scaffold_rl_translate_regions(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    matches_t* const matches) {
  // Paramters
  uint32_t* const rl_key_runs_acc = pattern->rl_runs_acc;
  uint32_t* const rl_text_runs_acc = text_trace->rl_runs_acc;
  const bool left_gap_alignment = (text_trace->strand==Forward);
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  // Translate all regions
  uint64_t i;
  for (i=0;i<num_alignment_regions;++i) {
    match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + i;
    // Translate into Text-Space
    match_alignment_region_set_type(match_alignment_region,match_alignment_region_approximate);
    // Translate CIGAR
    match_alignment_region_rl_translate(
        match_alignment_region,rl_key_runs_acc,rl_text_runs_acc,
        left_gap_alignment,matches->cigar_vector);
    // Region boundaries
    const uint64_t region_key_begin = match_alignment_region_get_key_begin(match_alignment_region);
    const uint64_t region_key_end = match_alignment_region_get_key_end(match_alignment_region);
    const uint64_t region_text_begin = match_alignment_region_get_text_begin(match_alignment_region);
    const uint64_t region_text_end = match_alignment_region_get_text_end(match_alignment_region);
    // Translate offsets
    match_alignment_region_set_key_begin(match_alignment_region,
        archive_text_rl_get_decoded_offset_exl(rl_key_runs_acc,region_key_begin));
    match_alignment_region_set_key_end(match_alignment_region,
        archive_text_rl_get_decoded_offset_exl(rl_key_runs_acc,region_key_end));
    match_alignment_region_set_text_begin(match_alignment_region,
        archive_text_rl_get_decoded_offset_exl(rl_text_runs_acc,region_text_begin));
    match_alignment_region_set_text_end(match_alignment_region,
        archive_text_rl_get_decoded_offset_exl(rl_text_runs_acc,region_text_end));
  }
}
/*
 * Region Chain Scaffolding
 */
void match_scaffold_region_chain(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    matches_t* const matches,
    mm_allocator_t* const mm_allocator) {
  // Translate alignment-regions
  const bool rl_space = match_scaffold->alignment_regions_rl;
  if (rl_space) {
    match_scaffold_rl_translate_regions(
        match_scaffold,pattern,text_trace,matches);
    match_scaffold->alignment_regions_rl = false;
  }
  // PROF
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_REGIONS,match_scaffold->num_alignment_regions);
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_COVERAGE,
      (100*match_scaffold->scaffolding_coverage)/pattern->key_length);
  // Scaffold chaining alignment-regions (from region-profile)
  match_scaffold_chain(match_scaffold,pattern,text_trace,!rl_space,mm_allocator);
  // PROF
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE,
      (100*match_scaffold->scaffolding_coverage)/pattern->key_length);
}
/*
 * Adaptive Scaffolding of the alignment (Best effort)
 */
void match_scaffold_adaptive(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    const uint64_t global_min_identity,
    const uint64_t matching_min_length,
    matches_t* const matches,
    mm_allocator_t* const mm_allocator) {
  PROFILE_START(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
  // Region-Chain scaffolding
  if (match_scaffold->scaffold_type == scaffold_none) {
    match_scaffold_region_chain(match_scaffold,pattern,text_trace,matches,mm_allocator);
    match_scaffold->scaffold_type = scaffold_region_chain; // Set scaffold approach
  }
  // Levenshtein scaffolding
  const uint64_t key_length = pattern->key_length;
  if (match_scaffold->scaffold_type <= scaffold_region_chain) {
    // Check coverage
    const uint64_t match_distance_min_bound = alignment->distance_min_bound;
    const uint64_t max_coverage_bound = BOUNDED_SUBTRACTION(key_length,match_distance_min_bound,0);
    if (max_coverage_bound < global_min_identity ||
        match_scaffold->scaffolding_coverage < max_coverage_bound) {
      // Scaffold from Levenshtein-alignment
      match_scaffold_levenshtein(
          match_scaffold,pattern,text_trace,
          alignment,matching_min_length,matches,mm_allocator);
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_COVERAGE,
          (100*match_scaffold->scaffolding_coverage)/key_length);
    }
    match_scaffold->scaffold_type = scaffold_levenshtein;
  }
  //    // DEBUG
  //    match_scaffold_print_pretty(stderr,matches,match_scaffold,
  //        align_input->key,align_input->key_length,
  //        align_input->text,align_input->text_length,mm_allocator);
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold (scaffold alignment)\n");
    tab_global_inc();
    match_scaffold_print(stderr,matches,match_scaffold);
    tab_global_dec();
  }
  // CHECK
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD_CHECK) {
    match_scaffold_check(match_scaffold,pattern,text_trace,matches);
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
}
/*
 * Check
 */
void match_scaffold_check(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    matches_t* const matches) {
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  uint64_t i;
  for (i=0;i<num_alignment_regions;++i) {
    match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + i;
    match_alignment_region_check(
        match_alignment_region,pattern->key,
        text_trace->text,matches->cigar_vector);
  }
}
/*
 * Sorting
 */
#define VECTOR_SORT_NAME                 match_alignment_region
#define VECTOR_SORT_TYPE                 match_alignment_region_t
#define VECTOR_SORT_CMP(a,b)             match_alignment_region_cmp_text_position(a,b)
#include "utils/vector_sort.h"
void match_scaffold_sort_alignment_regions(
    match_scaffold_t* const match_scaffold) {
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
    case scaffold_region_chain: tab_fprintf(stream,"  => Scaffold.type RegionChain\n"); break;
    case scaffold_levenshtein: tab_fprintf(stream,"  => Scaffold.type Levenshtein\n"); break;
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

