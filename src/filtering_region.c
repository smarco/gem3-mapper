/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_verify.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region.h"
#include "align.h"
//#include "align_bpm_distance.h"
//#include "match_scaffold.h"
//#include "match_align.h"
//#include "match_align_local.h"
//#include "match_align_dto.h"
//#include "output_map.h"

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
      "align-distance=%"PRId64","
      "matching-regions=%"PRIu64","
      "align-match=(%"PRIu64",%"PRIu64"))\n",
      filtering_region_status_label[region->status],region->begin_position,region->end_position,
      region->end_position-region->begin_position,
      region->align_distance==ALIGN_DISTANCE_INF ? (int64_t)-1 : (int64_t)region->align_distance,
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
