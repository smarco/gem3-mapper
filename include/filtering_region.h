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
#include "approximate_search_parameters.h"
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
  filtering_region_pending,               // On-hold state (unverified)
  filtering_region_unverified,            // Region still unverified
  filtering_region_accepted,              // Region accepted by pre-filters
  filtering_region_accepted_subdominant,  // Region accepted by pre-filters but sub-dominant
  filtering_region_aligned,               // Region aligned
  filtering_region_aligned_subdominant,   // Region aligned but sub-dominant in the end
  filtering_region_discarded              // Region discarded
} filtering_region_status_t;
typedef struct {
  // State
  filtering_region_status_t status;
  // Text-trace
  uint64_t text_trace_offset;
  // Location
  uint64_t begin_position;
  uint64_t effective_begin_position;
  uint64_t effective_end_position;
  // Regions Matching
  uint64_t num_regions_matching;
  region_matching_t* regions_matching;
  uint64_t coverage;
  // Alignment distance
  uint64_t align_distance;            // Distance (In case of approximate: max-bound)
  uint64_t align_distance_min_bound;  // Approximated distance (min-bound)
  uint64_t align_match_begin_column;  // Match begin column (inclusive)
  uint64_t align_match_end_column;    // Matching column (inclusive)
} filtering_region_t;
typedef struct {
  // State
  filtering_region_status_t status;
  // Location
  uint64_t effective_begin_position;
  uint64_t effective_end_position;
} verified_region_t;

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
GEM_INLINE bool filtering_region_align(
    filtering_region_t* const filtering_region,const text_collection_t* const candidates_collection,
    search_parameters_t* const search_parameters,const strand_t search_strand,
    const pattern_t* const pattern,matches_t* const matches,match_trace_t* const match_trace,
    mm_stack_t* const mm_stack);

/*
 * Verify
 */
GEM_INLINE bool filtering_region_verify(
    filtering_region_t* const filtering_region,
    const text_collection_t* const candidates_collection,
    search_parameters_t* const search_parameters,const pattern_t* const pattern);
GEM_INLINE uint64_t filtering_region_verify_extension(
    vector_t* const filtering_regions,
    const text_collection_t* const candidates_collection,
    const uint64_t const text_trace_offset,const uint64_t index_position,
    search_parameters_t* const search_parameters,const pattern_t* const pattern);

/*
 * Display
 */
GEM_INLINE void filtering_region_print_matching_regions(
    FILE* const stream,region_matching_t* const regions_matching,const uint64_t num_regions_matching,
    const uint64_t begin_position,const uint64_t end_position,const uint64_t filtering_region_idx);

#endif /* FILTERING_REGION_H_ */
