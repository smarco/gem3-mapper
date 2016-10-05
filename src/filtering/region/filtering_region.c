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
 * DESCRIPTION:
 *   Filtering region data structure holds all the required information
 *   of a candidate region of the index that can potentially align against
 *   the search pattern
 */

#include "align/alignment.h"
#include "text/pattern.h"
#include "filtering/region/filtering_region.h"
#include "archive/archive_text_rl.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION                GEM_DEEP_DEBUG
#define DEBUG_DISPLAY_MATCHING_TEXT

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
 * Retrieve filtering region text-candidate
 */
void filtering_region_retrieve_text(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    archive_text_t* const archive_text,
    text_collection_t* const text_collection) {
  // Check already retrieved
  if (filtering_region->text_trace_offset != UINT64_MAX) return;
  // Retrieve Text
  uint64_t text_length;
  if (archive_text->run_length) {
    // Translate RL-Text positions (RL-text encoded)
    const uint64_t text_begin_position = archive_text_rl_position_translate(
        archive_text,filtering_region->text_begin_position,text_collection->mm_text);
    const uint64_t text_end_position = archive_text_rl_position_translate(
        archive_text,filtering_region->text_end_position,text_collection->mm_text);
    // Retrieve Text
    text_length = text_end_position - text_begin_position;
    const uint64_t text_trace_offset = archive_text_retrieve_collection(
        archive_text,text_collection,text_begin_position,text_length,false,true);
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
    filtering_region->text_trace_offset = archive_text_retrieve_collection(
        archive_text,text_collection,text_position,text_length,false,false);
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
      pattern_trimmed_init(pattern,
          &filtering_region->bpm_pattern_trimmed,
          &filtering_region->bpm_pattern_trimmed_tiles,
          filtering_region->key_trimmed_length,
          filtering_region->key_trim_left,mm_stack);
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
void filtering_region_print_region_text(
    FILE* const stream,
    filtering_region_t* const region,
    const text_collection_t* const text_collection) {
  if (text_collection!=NULL && region->text_trace_offset!=UINT64_MAX) {
    // Retrieve text
    const uint64_t text_length = region->text_end_position-region->text_begin_position;
    const text_trace_t* const text_trace = text_collection_get_trace(text_collection,region->text_trace_offset);
    uint8_t* const text = text_trace->text;
    // Allocate display text
    const uint64_t max_printed_length = MIN(200,text_length);
    uint64_t i;
#ifdef DEBUG_DISPLAY_MATCHING_TEXT
    char* const display_text = malloc(max_printed_length);
    uint64_t s, p;
    for (i=0;i<max_printed_length;++i) display_text[i] = 'a'+(dna_decode(text[i])-'A');
    // Upper-case alignment-regions
    match_scaffold_t* const match_scaffold = &region->match_scaffold;
    for (s=0;s<match_scaffold->num_alignment_regions;++s) {
      match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + s;
      const uint64_t region_text_end = match_alignment_region_get_text_end(match_alignment_region);
      const uint64_t region_text_begin = match_alignment_region_get_text_begin(match_alignment_region);
      const uint64_t max_text_scope = MIN(max_printed_length,region_text_end);
      for (p=region_text_begin;p<max_text_scope;++p) display_text[p] = dna_decode(text[p]);
    }
    // Display
    tab_fprintf(stream,"  => Text %.*s\n",max_printed_length,display_text);
    // Free
    free(display_text);
#else
    tab_fprintf(stream,"  => Text ");
    for (i=0;i<max_printed_length;++i) {
      fprintf(stream,"%c",dna_decode(text[i]));
    }
    fprintf(stream,"\n");
#endif
  }
}
void filtering_region_print_alignment(
    FILE* const stream,
    alignment_t* const alignment) {
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  if (alignment_tiles==NULL) return;
  uint64_t i;
  tab_fprintf(stream,"  => Alignment-Tiles total=%lu dist-bound=%lu",
      alignment->num_tiles,
      alignment->distance_min_bound);
  for (i=0;i<alignment->num_tiles;++i) {
    fprintf(stream," (%lu,%lu][%lu]",
        alignment_tiles[i].text_begin_offset,
        alignment_tiles[i].text_end_offset,
        alignment_tiles[i].distance);
  }
  fprintf(stream,"\n");
}
void filtering_region_print(
    FILE* const stream,
    filtering_region_t* const region,
    const text_collection_t* const text_collection,
    const bool print_region_text,
    const bool print_alignment_regions,
    const bool print_alignment) {
  alignment_t* const alignment = &region->alignment;
  tab_fprintf(stream,"  => Region %s [%"PRIu64",%"PRIu64") "
      "(total-bases=%"PRIu64","
      "alignment-regions=%"PRIu64","
      "align-distance=%"PRId64",",
      filtering_region_status_label[region->status],
      region->text_begin_position,region->text_end_position,
      region->text_end_position-region->text_begin_position,
      region->match_scaffold.num_alignment_regions,
      alignment->distance_min_bound==ALIGN_DISTANCE_INF ?
          (int64_t)-1 : (int64_t)alignment->distance_min_bound);
  if (alignment->distance_min_bound!=ALIGN_DISTANCE_INF) {
    alignment_tile_t* const alignment_tile = alignment->alignment_tiles;
    fprintf(stream,"align-range=(%"PRIu64",%"PRIu64"))\n",
        alignment_tile[0].text_begin_offset,
        alignment_tile[alignment->num_tiles-1].text_end_offset);
  } else {
    fprintf(stream,"align-range=n/a)\n");
  }
  // Print Region-Text
  if (print_region_text) {
    if (text_collection!=NULL && region->text_trace_offset!=UINT64_MAX) {
      tab_global_inc();
      filtering_region_print_region_text(stream,region,text_collection);
      tab_global_dec();
    }
  }
  // Print Alignment-Regions
  if (print_alignment_regions) {
    tab_global_inc();
    match_scaffold_print(stream,NULL,&region->match_scaffold);
    tab_global_dec();
  }
  // Alignment
  if (print_alignment) {
    tab_global_inc();
    filtering_region_print_alignment(stream,alignment);
    tab_global_dec();
  }
}
