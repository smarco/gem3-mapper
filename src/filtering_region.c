/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_verify.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region.h"
#include "align.h"
#include "pattern.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION  GEM_DEEP_DEBUG
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
    [6] = "aligned-subdominant",
    [7] = "aligned-unbounded"
};

/*
 * Accessors
 */
void filtering_region_add(
    vector_t* const filtering_regions,const uint64_t text_trace_offset,
    const uint64_t begin_position,const uint64_t end_position,
    const uint64_t align_distance,const uint64_t align_match_end_column) {
  filtering_region_t* filtering_region;
  vector_alloc_new(filtering_regions, filtering_region_t, filtering_region);
  // State
  filtering_region->status = filtering_region_accepted;
  // Text-trace
  filtering_region->text_trace_offset = text_trace_offset;
  // Location
  filtering_region->begin_position = begin_position;
  filtering_region->end_position = end_position;
  filtering_region->base_begin_position_offset = 0;
  filtering_region->base_end_position_offset = end_position-begin_position;
  filtering_region->key_trim_left = 0;
  filtering_region->key_trim_right = 0;
  // Trimmed Pattern
  filtering_region->bpm_pattern_trimmed = NULL;
  filtering_region->bpm_pattern_trimmed_tiles = NULL;
  // Regions Matching
  match_scaffold_init(&filtering_region->match_scaffold);
  // Alignment distance
  filtering_region->align_distance = align_distance;
  filtering_region->align_match_end_column = align_match_end_column;
}
/*
 * Translate filtering-region offsets to RL-space
 */
void filtering_region_translate_run_length(
    filtering_region_t* const filtering_region,
    text_collection_t* const text_collection,pattern_t* const pattern) {
//  // Retrieve text
//  text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
//  // Keep original 'effective begin position'
//  filtering_region->source_begin_position = filtering_region->begin_position;
//  // Translate key-offsets
//  uint8_t* const rl_key_runs = pattern->rl_runs;
//  const uint64_t rl_key_length = pattern->rl_key_length;
//  filtering_region->key_trim_left = archive_text_rl_encode(rl_key_runs,rl_key_length,filtering_region->key_trim_left);
//  filtering_region->key_trim_right = archive_text_rl_encode(rl_key_runs,rl_key_length,filtering_region->key_trim_right);

//  // Translate text-positions
//  uint8_t* const rl_text_runs = text_trace->rl_runs;
//  const uint64_t rl_text_length = text_trace->rl_text_length;
//  filtering_region->begin_position = archive_text_rl_encode(rl_text_runs,rl_text_length,filtering_region->begin_position);
//  filtering_region->end_position = archive_text_rl_encode(rl_text_runs,rl_text_length,filtering_region->end_position);
//  // Translate text-offsets
//  const uint64_t base_begin_position =
//      archive_text_rl_encode(rl_text_runs,rl_text_length,filtering_region->base_begin_position_offset);

//  filtering_region->base_begin_position_offset = ;
//  filtering_region->base_end_position_offset = archive_text_rl_encode(rl_text_runs,rl_text_length,filtering_region->base_end_position_offset);;
}
/*
 * Retrieve filtering region text-candidate
 */
void filtering_region_retrieve_text(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    text_collection_t* const text_collection,mm_stack_t* const mm_stack) {
  // Check already retrieved
  if (filtering_region->text_trace_offset != UINT64_MAX) return;
  // Retrieve Text
  const uint64_t text_position = filtering_region->begin_position;
  const uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
  const bool run_length_text = archive_text->run_length;
  filtering_region->text_trace_offset =
      archive_text_retrieve_collection(archive_text,text_collection,
          text_position,text_length,false,run_length_text,mm_stack);
  // Adapt RL-Text
  if (run_length_text) {

  }
}
/*
 * Sorting
 */
int filtering_region_locator_cmp_position(
    const filtering_region_locator_t* const a,const filtering_region_locator_t* const b) {
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
    FILE* const stream,filtering_region_t* const region,
    const text_collection_t* const text_collection,const bool print_matching_regions) {
  tab_fprintf(stream,"  => Region %s [%"PRIu64",%"PRIu64") "
      "(total-bases=%"PRIu64","
      "align-distance=(%"PRId64",%"PRId64"),"
      "matching-regions=%"PRIu64","
      "align-match=(%"PRIu64",%"PRIu64"))\n",
      filtering_region_status_label[region->status],region->begin_position,region->end_position,
      region->end_position-region->begin_position,
      region->align_distance==ALIGN_DISTANCE_INF ? (int64_t)-1 : (int64_t)region->align_distance,
      region->align_distance_min_bound==ALIGN_DISTANCE_INF ? (int64_t)-1 : (int64_t)region->align_distance_min_bound,
      region->match_scaffold.num_scaffold_regions,
      region->align_match_begin_column,region->align_match_end_column);
  if (text_collection!=NULL) {
    if (region->text_trace_offset == UINT64_MAX) {
      tab_fprintf(stream,"    => Text 'n/a'\n");
    } else {
      // Retrieve text
      const uint64_t text_length = region->end_position-region->begin_position;
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
  }
  if (print_matching_regions) {
    tab_global_inc();
    match_scaffold_print(stream,NULL,&region->match_scaffold);
    tab_global_dec();
  }
}
