/*
 * PROJECT: GEMMapper
 * FILE: filtering_region.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_H_
#define FILTERING_REGION_H_

#include "essentials.h"

#include "text_collection.h"
#include "matches_align.h"
#include "pattern.h"
#include "bpm_align.h"

/*
 * Debug
 */
#define DEBUG_ALIGN_CANDIDATES  false

/*
 * Filtering Region
 */
typedef enum {
  filtering_region_new,
  filtering_region_verified,
  filtering_region_accepted,
  filtering_region_accepted_subdominant,
  filtering_region_aligned,
  filtering_region_discarded
} filtering_region_status_t;
typedef struct {
  // State
  filtering_region_status_t status;
  // Text-trace
  uint64_t text_trace_offset;
  // Location
  uint64_t begin_position;
  uint64_t effective_end_position;
  uint64_t effective_begin_position;
  // Regions Matching
  uint64_t num_regions_matching;
  region_matching_t* regions_matching;
  uint64_t coverage;
  // Alignment distance
  uint64_t align_distance;
  uint64_t align_match_column;
} filtering_region_t;

/*
 * Accessors
 */
GEM_INLINE text_trace_t* filtering_region_get_text_trace(
    const filtering_region_t* const filtering_region,const text_collection_t* const candidates_collection);

/*
 * Sorting
 */
GEM_INLINE void filtering_region_sort_regions_matching(const filtering_region_t* const filtering_region);

/*
 * Matching regions
 */
GEM_INLINE void filtering_region_exact_extend(
    filtering_region_t* const filtering_region,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const bool* const allowed_enc);
GEM_INLINE void filtering_region_chain_matching_regions(
    filtering_region_t* const filtering_region,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const bool* const allowed_enc,
    const uint64_t max_error,mm_stack_t* const stack);

/*
 * (Re)Align
 */
GEM_INLINE void filtering_region_align(
    filtering_region_t* const accepted_region,const text_collection_t* const candidates_collection,
    const alignment_model_t alignment_model,const bool* const allowed_enc,
    const swg_penalties_t* const swg_penalties,const strand_t search_strand,
    const pattern_t* const pattern,const uint8_t* const key,const uint64_t key_length,
    matches_t* const matches,mm_stack_t* const mm_stack);

/*
 * Verify
 */
GEM_INLINE bool filtering_region_verify(
    filtering_region_t* const candidate_region,const text_collection_t* const candidates_collection,
    const alignment_model_t alignment_model,const bool* const allowed_enc,const pattern_t* const pattern,
    const uint8_t* const key,const uint64_t key_length,const uint64_t max_effective_filtering_error);

/*
 * Display
 */
GEM_INLINE void filtering_region_print_matching_regions(
    FILE* const stream,filtering_region_t* const filtering_region,const uint64_t filtering_region_idx);

#endif /* FILTERING_REGION_H_ */
