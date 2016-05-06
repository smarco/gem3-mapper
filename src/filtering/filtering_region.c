/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_verify.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_region.h"
#include "align/align.h"
#include "data_structures/pattern.h"
#include "archive/archive_text_rl.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION                GEM_DEEP_DEBUG
#define DEBUG_FILTERING_REGION_DISPLAY_TEXT_MATCHING_REGIONS

/*
 * Data Types & Tables
 */
const char* filtering_region_status_label[] =
{
    [0] = "pending",
    [1] = "unverified",
    [2] = "verified-discarded",
    [3] = "accepted",
    [4] = "accepted-subdominant",
    [5] = "aligned",
};
/*
 * Accessors
 */
void filtering_region_add(
    vector_t* const filtering_regions,
    const uint64_t text_trace_offset,
    const uint64_t begin_position,
    const uint64_t end_position,
    const uint64_t align_distance,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset,
    mm_stack_t* const mm_stack) {
  filtering_region_t* filtering_region;
  vector_alloc_new(filtering_regions, filtering_region_t, filtering_region);
  // State
  filtering_region->status = filtering_region_accepted;
  // Text-trace
  filtering_region->text_trace_offset = text_trace_offset;
  // Location
  filtering_region->text_begin_position = begin_position;
  filtering_region->text_end_position = end_position;
  filtering_region->key_trim_left = 0;
  filtering_region->key_trim_right = 0;
  // Trimmed Pattern
  filtering_region->bpm_pattern_trimmed = NULL;
  filtering_region->bpm_pattern_trimmed_tiles = NULL;
  // Regions Matching
  match_scaffold_init(&filtering_region->match_scaffold);
  // Alignment distance
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  region_alignment->distance_min_bound = align_distance;
  region_alignment->num_tiles = 1;
  region_alignment_tile_t* const alignment_tiles = mm_stack_calloc(mm_stack,1,region_alignment_tile_t,false);
  region_alignment->alignment_tiles = alignment_tiles;
  alignment_tiles->match_distance = align_distance;
  alignment_tiles->text_begin_offset = text_begin_offset;
  alignment_tiles->text_end_offset = text_end_offset;
}
/*
 * Retrieve filtering region text-candidate
 */
void filtering_region_retrieve_text(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    archive_text_t* const archive_text,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack) {
  // Check already retrieved
  if (filtering_region->text_trace_offset != UINT64_MAX) return;
  // Retrieve Text
  uint64_t text_length;
  if (archive_text->run_length) {
    // Translate RL-Text positions (RL-text encoded)
    const uint64_t text_begin_position = archive_text_rl_position_translate(
        archive_text,filtering_region->text_begin_position,mm_stack);
    const uint64_t text_end_position = archive_text_rl_position_translate(
        archive_text,filtering_region->text_end_position,mm_stack);
    // Retrieve Text
    text_length = text_end_position - text_begin_position;
    const uint64_t text_trace_offset = archive_text_retrieve_collection(
        archive_text,text_collection,text_begin_position,text_length,false,true,mm_stack);
    // Configure filtering-region
    text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
    filtering_region->text_trace_offset = text_trace_offset;
    filtering_region->text_begin_position = text_begin_position;
    filtering_region->text_end_position = text_end_position;
    // Set fix position
    filtering_region->text_source_region_offset = archive_text_rl_get_decoded_offset_exl(
        text_trace->rl_runs_acc,filtering_region->text_source_region_offset);
    filtering_region->key_source_region_offset = archive_text_rl_get_decoded_offset_exl(
        pattern->rl_runs_acc,filtering_region->key_source_region_offset);
    // Compute key trims
    filtering_region_compute_key_trims(filtering_region,pattern);
    //    // DEBUG
    //    text_trace_t* const text_trace = text_collection_get_trace(
    //        text_collection,filtering_region->text_trace_offset);
    //    fprintf(stderr,">Text\n");
    //    dna_buffer_print(stderr,text_trace->text,text_trace->text_length,false);
    //    fprintf(stderr,"\n");
    //    fprintf(stderr,">RL.Text\n");
    //    dna_buffer_print(stderr,text_trace->rl_text,text_trace->rl_text_length,false);
    //    fprintf(stderr,"\n");
    //    uint64_t i; for (i=0;i<text_trace->rl_text_length;++i) fprintf(stderr,"[%d]",text_trace->rl_runs[i]);
    //    fprintf(stderr,"\n");
  } else {
    // Retrieve Text
    const uint64_t text_position = filtering_region->text_begin_position;
    text_length = filtering_region->text_end_position - filtering_region->text_begin_position;
    filtering_region->text_trace_offset =
        archive_text_retrieve_collection(archive_text,text_collection,
            text_position,text_length,false,false,mm_stack);
  }
}
/*
 * Prepare Alignment
 */
void filtering_region_alignment_prepare(
    filtering_region_t* const filtering_region,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    mm_stack_t* const mm_stack) {
  // Check region-alignment
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  if (region_alignment->alignment_tiles!=NULL) return; // Already initialized
  // Allocate region-alignment
  const uint64_t num_tiles = bpm_pattern_tiles->num_pattern_tiles;
  region_alignment->num_tiles = num_tiles;
  region_alignment->distance_min_bound = bpm_pattern->pattern_length;
  region_alignment->alignment_tiles = mm_stack_calloc(mm_stack,num_tiles,region_alignment_tile_t,false);
  // Init all tiles
  const uint64_t text_length = filtering_region->text_end_position-filtering_region->text_begin_position;
  region_alignment_tile_t* const alignment_tiles = region_alignment->alignment_tiles;
  if (num_tiles==1) {
    alignment_tiles->match_distance = ALIGN_DISTANCE_INF;
    alignment_tiles->text_end_offset = text_length;
    alignment_tiles->text_begin_offset = 0;
  } else {
    // Calculate tile dimensions
    const uint64_t max_error = MIN(filtering_region->max_error,bpm_pattern->pattern_length);
    pattern_tiled_t pattern_tiled;
    pattern_tiled_init(&pattern_tiled,bpm_pattern->pattern_length,
        bpm_pattern_tiles->tile_length,text_length,max_error);
    uint64_t tile_pos;
    for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
      // Init Tile
      alignment_tiles[tile_pos].match_distance = ALIGN_DISTANCE_INF;
      alignment_tiles[tile_pos].text_end_offset = pattern_tiled.tile_offset+pattern_tiled.tile_wide;
      alignment_tiles[tile_pos].text_begin_offset = pattern_tiled.tile_offset;
      // Calculate next tile
      pattern_tiled_calculate_next(&pattern_tiled);
    }
  }
}
/*
 * Compute key trims
 */
void filtering_region_compute_key_trims(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  // Compute key trims
  const uint64_t key_length = pattern->key_length;
  const uint64_t text_fix_begin = filtering_region->text_source_region_offset;
  const uint64_t key_fix_begin = filtering_region->key_source_region_offset;
  const uint64_t text_length = filtering_region->text_end_position - filtering_region->text_begin_position;
  if (pattern->key_length > text_length) {
    // Compute trim offsets
    filtering_region->key_trim_left = (key_fix_begin > text_fix_begin) ? key_fix_begin - text_fix_begin : 0;
    const uint64_t text_fix_end = text_length - text_fix_begin;
    const uint64_t key_fix_end = key_length - key_fix_begin;
    filtering_region->key_trim_right = (key_fix_end > text_fix_end) ? key_fix_end - text_fix_end : 0;
    filtering_region->key_trimmed_length =
        pattern->key_length - filtering_region->key_trim_left - filtering_region->key_trim_right;
    // Set max-error
    const double max_error_factor = (double)pattern->max_effective_filtering_error / (double)key_length;
    const double max_bandwidth_factor = (double)pattern->max_effective_bandwidth / (double)key_length;
    filtering_region->max_error = max_error_factor * (double)filtering_region->key_trimmed_length;
    filtering_region->max_bandwidth = max_bandwidth_factor * (double)filtering_region->key_trimmed_length;
    // Set trimmed & init fields
    filtering_region->key_trimmed = true;
    filtering_region->bpm_pattern_trimmed = NULL;
    filtering_region->bpm_pattern_trimmed_tiles = NULL;
  } else {
    filtering_region->key_trim_left = 0;
    filtering_region->key_trim_right = 0;
    filtering_region->key_trimmed = false;
  }
}
/*
 * Select proper BPM-Pattern
 */
void filtering_region_bpm_pattern_select(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    bpm_pattern_t** const bpm_pattern,
    bpm_pattern_t** const bpm_pattern_tiles,
    mm_stack_t* const mm_stack) {
  // Select BPM-Pattern
  if (filtering_region->key_trimmed) {
    // Check compiled
    if (filtering_region->bpm_pattern_trimmed==NULL) {
      pattern_trimmed_init(pattern,&filtering_region->bpm_pattern_trimmed,
          &filtering_region->bpm_pattern_trimmed_tiles,filtering_region->key_trimmed_length,
          filtering_region->key_trim_left,filtering_region->key_trim_right,mm_stack);
    }
    *bpm_pattern = filtering_region->bpm_pattern_trimmed;
    *bpm_pattern_tiles = filtering_region->bpm_pattern_trimmed_tiles;
  } else {
    *bpm_pattern = pattern->bpm_pattern;
    *bpm_pattern_tiles = pattern->bpm_pattern_tiles;
  }
}
/*
 * Sorting
 */
int filtering_region_locator_cmp_position(
    const filtering_region_locator_t* const a,
    const filtering_region_locator_t* const b) {
  return a->position - b->position;
}
void filtering_region_locator_sort_positions(vector_t* const filtering_region_locators) {
  void* array = vector_get_mem(filtering_region_locators,filtering_region_locator_t);
  const size_t count = vector_get_used(filtering_region_locators);
  qsort(array,count,sizeof(filtering_region_locator_t),
      (int (*)(const void *,const void *))filtering_region_locator_cmp_position);
}
/*
 * Display
 */
void filtering_region_print(
    FILE* const stream,
    filtering_region_t* const region,
    const text_collection_t* const text_collection,
    const bool print_matching_regions) {
  region_alignment_t* const region_alignment = &region->region_alignment;
  tab_fprintf(stream,"  => Region %s [%"PRIu64",%"PRIu64") "
      "(total-bases=%"PRIu64","
      "scaffold-regions=%"PRIu64","
      "align-distance=%"PRId64",",
      filtering_region_status_label[region->status],
      region->text_begin_position,region->text_end_position,
      region->text_end_position-region->text_begin_position,
      region->match_scaffold.num_scaffold_regions,
      region_alignment->distance_min_bound==ALIGN_DISTANCE_INF ?
          (int64_t)-1 : (int64_t)region_alignment->distance_min_bound);
  if (region_alignment->distance_min_bound!=ALIGN_DISTANCE_INF) {
    region_alignment_tile_t* const alignment_tile = region_alignment->alignment_tiles;
    fprintf(stream,"align-range=(%"PRIu64",%"PRIu64"))\n",
        alignment_tile[0].text_begin_offset,
        alignment_tile[region_alignment->num_tiles-1].text_end_offset);
  } else {
    fprintf(stream,"align-range=n/a)\n");
  }
  if (text_collection!=NULL && region->text_trace_offset != UINT64_MAX) {
    // Retrieve text
    const uint64_t text_length = region->text_end_position-region->text_begin_position;
    const text_trace_t* const text_trace = text_collection_get_trace(text_collection,region->text_trace_offset);
    uint8_t* const text = text_trace->text;
    // Allocate display text
    const uint64_t max_printed_length = MIN(200,text_length);
    uint64_t i;
#ifdef DEBUG_FILTERING_REGION_DISPLAY_TEXT_MATCHING_REGIONS
    char* const display_text = malloc(max_printed_length);
    uint64_t s, p;
    for (i=0;i<max_printed_length;++i) display_text[i] = 'a'+(dna_decode(text[i])-'A');
    // Upper-case matching regions
    match_scaffold_t* const match_scaffold = &region->match_scaffold;
    for (s=0;s<match_scaffold->num_scaffold_regions;++s) {
      region_matching_t* const region_matching = match_scaffold->scaffold_regions + s;
      const uint64_t max_text_scope = MIN(max_printed_length,region_matching->text_end);
      for (p=region_matching->text_begin;p<max_text_scope;++p) display_text[p] = dna_decode(text[p]);
    }
    // Display
    tab_fprintf(stream,"    => Text %.*s\n",max_printed_length,display_text);
    // Free
    free(display_text);
#else
    tab_fprintf(stream,"    => Text ");
    for (i=0;i<max_printed_length;++i) {
      fprintf(stream,"%c",dna_decode(text[i]));
    }
    fprintf(stream,"\n");
#endif
  }
  if (print_matching_regions) {
    tab_global_inc();
    match_scaffold_print(stream,NULL,&region->match_scaffold);
    tab_global_dec();
  }
}
