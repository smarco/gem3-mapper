/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_align.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates_align_local.h"
#include "filtering/filtering_candidates_align.h"
#include "filtering/filtering_region_verify.h"
#include "align/align.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Exclude tiles
 */
void filtering_candidates_align_local_exclude_tiles(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    mm_stack_t* const mm_stack) {
  // Parameters
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  // Get BPM-Pattern
  bpm_pattern_t* bpm_pattern, *bpm_pattern_tiles;
  filtering_region_bpm_pattern_select(filtering_region,
      pattern,&bpm_pattern,&bpm_pattern_tiles,mm_stack);
  // Prepare Alignment
  filtering_region_alignment_prepare(filtering_region,
      bpm_pattern,bpm_pattern_tiles,mm_stack);
  // Check number of tiles
  const uint64_t num_pattern_tiles = bpm_pattern_tiles->num_pattern_tiles;
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  region_alignment_tile_t* const alignment_tiles = region_alignment->alignment_tiles;
  if (num_pattern_tiles > 1) {
    // Sort matching regions by text-offsets
    match_scaffold_sort_regions_matching(match_scaffold);
    // Exclude tiles without supporting matching-region
    region_matching_t* region_matching = match_scaffold->scaffold_regions;
    uint64_t tile_pos, tile_key_begin = 0, region_pos;
    for (tile_pos=0,region_pos=0;tile_pos<num_pattern_tiles;++tile_pos) {
      bpm_pattern_t* const bpm_pattern_tile = bpm_pattern_tiles+tile_pos;
      const uint64_t tile_key_end = tile_key_begin + bpm_pattern_tile->pattern_length;
      // Skip resolved tiles
      if (alignment_tiles[tile_pos].match_distance!=ALIGN_DISTANCE_INF) continue;
      // Find nearest matching region
      while (region_matching!=NULL && region_matching->key_end <= tile_key_begin) {
        region_matching = (++region_pos < num_scaffold_regions) ? region_matching+1 : NULL;
      }
      // Check overlapping
      if (region_matching==NULL || tile_key_end <= region_matching->key_begin) {
        alignment_tiles[tile_pos].match_distance = ALIGN_DISABLED;
      }
      // Next
      tile_key_begin = tile_key_end;
    }
  }
}
/*
 * Filtering Candidates (Re)Alignment
 */
void filtering_candidates_align_local(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    matches_t* const matches) {
  PROFILE_START(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Add pending local matches (found so far)
  locator_t* const locator = filtering_candidates->archive->locator;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  matches_add_pending_local_matches(matches,locator,mm_stack);
  // Check total alignments found
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  const uint64_t max_reported_matches = select_parameters->max_reported_matches;
  uint64_t total_matches = matches_get_num_match_traces(matches);
  if (total_matches < max_reported_matches) {
    // Clear cache
    filtering_region_cache_clear(&filtering_candidates->filtering_region_cache);
    // Sort by scaffold-coverage
    filtering_regions_sort_scaffold_coverage(filtering_candidates->filtering_regions);
    // Traverse all discarded regions & local-align
    const uint64_t num_regions = vector_get_used(filtering_candidates->discarded_regions);
    filtering_region_t* filtering_region = vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t);
    uint64_t i; // last_region_coverage;
    PROF_ADD_COUNTER(GP_CANDIDATE_REGION_LOCAL,num_regions);
    for (i=0;i<num_regions;++i,++filtering_region) {
//      // Check max-reported matches & coverage
//      const uint64_t scaffolding_coverage = filtering_region->match_scaffold.scaffolding_coverage;
//      if (total_matches >= max_reported_matches && last_region_coverage > scaffolding_coverage) break;
//      last_region_coverage = scaffolding_coverage;
      // Retrieve Text
      filtering_region_retrieve_text(
          filtering_region,pattern,filtering_candidates->archive->text,
          filtering_candidates->text_collection,mm_stack);
      // Exclude not-supported regions for local-alignment
      filtering_candidates_align_local_exclude_tiles(filtering_region,pattern,mm_stack);
      // Align Region
      PROF_INC_COUNTER(GP_CANDIDATE_REGION_LOCAL_ALIGNED);
      const bool match_added = filtering_candidates_align_region(
          filtering_candidates,filtering_region,pattern,emulated_rc_search,true,false,matches);
      if (match_added) ++total_matches;
    }
    // Update Used
    vector_clear(filtering_candidates->discarded_regions);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (local_align)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
}
