/*
 * PROJECT: GEMMapper
 * FILE: filtering_region.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_H_
#define FILTERING_REGION_H_

#include "utils/essentials.h"
#include "data_structures/sequence.h"
#include "data_structures/pattern.h"
#include "data_structures/text_collection.h"
#include "archive/archive_text.h"
#include "matches/match_scaffold.h"

/*
 * Filtering Region Status
 */
typedef enum {
  filtering_region_pending=0,               // On-hold state (unverified)
  filtering_region_unverified=1,            // Region still unverified
  filtering_region_verified_discarded=2,    // Region verified but too distant
  filtering_region_accepted=3,              // Region accepted by pre-filters
  filtering_region_accepted_subdominant=4,  // Region accepted by pre-filters but sub-dominant (by ranking)
  filtering_region_aligned=5,               // Region aligned
} filtering_region_status_t;
extern const char* filtering_region_status_label[6];
/*
 * Filtering Region
 */
typedef struct {
  /* State */
  filtering_region_status_t status;         // Filtering Region Status
  /* Source Region Offset */
  uint64_t text_source_region_offset;       // Text-Offset to the begin of the source-region
  uint64_t key_source_region_offset;        // Key-Offset to the begin of the source-region
  /* Text */
  uint64_t text_trace_offset;               // Text-trace
  uint64_t text_begin_position;             // Region effective begin position (adjusted to error boundaries)
  uint64_t text_end_position;               // Region effective end position (adjusted to error boundaries)
  /* Key */
  uint64_t key_trimmed_length;              // Key trimmed length
  uint64_t key_trim_left;                   // Key left trim  (due to position correction wrt the reference limits)
  uint64_t key_trim_right;                  // Key right trim (due to position correction wrt the reference limits)
  bool key_trimmed;                         // Key has been trimmed (due to dangling ends; projected to the text-candidate)
  bpm_pattern_t* bpm_pattern_trimmed;       // BPM-Pattern Trimmed
  bpm_pattern_t* bpm_pattern_trimmed_tiles; // BPM-Pattern Trimmed & Tiled
  /* Alignment */
  uint64_t max_error;                       // Maximum effective error
  uint64_t max_bandwidth;                   // Maximum effective bandwidth
  region_alignment_t region_alignment;      // Filtering Region Alignment
  match_scaffold_t match_scaffold;          // Matching regions supporting the filtering region (Scaffolding)
} filtering_region_t;
typedef struct {
  /* Source Region Offset */
  uint32_t text_source_region_offset;       // Text-Offset to the begin of the source-region
  uint32_t key_source_region_offset;        // Key-Offset to the begin of the source-region
  /* Text */
  uint64_t text_begin_position;             // Region effective begin position (adjusted to error boundaries)
  uint64_t text_end_position;               // Region effective end position (adjusted to error boundaries)
  /* Alignment */
  region_alignment_t region_alignment;      // Filtering Region Alignment
  region_matching_t* scaffold_regions;
  uint64_t num_scaffold_regions;
} filtering_region_buffered_t;
/*
 * Filtering Region Auxiliary Structures
 */
typedef struct {
  /* Filtering region */
  filtering_region_t* filtering_region;
  /* Location */
  sequence_end_t sequence_end;
  strand_t strand;
  uint64_t position;
} filtering_region_locator_t;
typedef struct {
  /* Location */
  uint64_t begin_position;
  uint64_t end_position;
} verified_region_t;

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
    mm_stack_t* const mm_stack);

/*
 * Retrieve filtering region text-candidate
 */
void filtering_region_retrieve_text(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    archive_text_t* const archive_text,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack);

/*
 * Prepare Alignment
 */
void filtering_region_alignment_prepare(
    filtering_region_t* const filtering_region,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    mm_stack_t* const mm_stack);

/*
 * Filtering Region Key Trims
 */
void filtering_region_compute_key_trims(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern);

/*
 * Filtering Region BPM-Pattern
 */
void filtering_region_bpm_pattern_select(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    bpm_pattern_t** const bpm_pattern,
    bpm_pattern_t** const bpm_pattern_tiles,
    mm_stack_t* const mm_stack);

/*
 * Sorting
 */
void filtering_region_sort_regions_matching(const filtering_region_t* const filtering_region);
void filtering_region_locator_sort_positions(vector_t* const filtering_region_locators);

/*
 * Display
 */
void filtering_region_print(
    FILE* const stream,
    filtering_region_t* const region,
    const text_collection_t* const text_collection,
    const bool print_matching_regions);

#endif /* FILTERING_REGION_H_ */
