/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_CANDIDATES_H_
#define FILTERING_CANDIDATES_H_

#include "essentials.h"
#include "interval_set.h"
#include "matches.h"

/*
 * Candidates
 */
typedef struct {
  // Ranges of the region [end,start)
  uint64_t start;
  uint64_t end;
  // Number of differences in the region
  uint64_t errors;
  // Internals (Overloaded fields)
  uint64_t overload_field_0;
  uint64_t overload_field_1;
  uint64_t overload_field_2;
} candidate_region_t;
typedef struct {
  // Region Info
  candidate_region_t* candidate_region;
  // Position
  uint64_t trace_offset;         // Trace-offset
  uint64_t begin_index_position; // Index position
  uint64_t distance;             // Levenshtein Score
} candidate_position_t;
/*
 * Filtering Candidates Vector
 */
typedef struct {
  // Candidates {Regions} pending to be verified
  svector_t* pending_candidate_regions;                    // (candidate_region_t)
  svector_iterator_t pending_candidate_regions_iterator;   // Writing Iterator (always appending)
  // Candidates Positions
  svector_t* pending_candidate_positions;                  // (candidate_position_t)
  svector_iterator_t pending_candidate_positions_iterator; // Writing Iterator
//  // Candidates Text-Collection (Stores candidates Text-blocks)
//  text_collection_t* candidates_collection;
  // Checked Positions
  svector_t* checked_candidate_positions;                  // (candidate_position_t)
  svector_iterator_t checked_candidate_positions_iterator; // Writing Iterator
} filtering_candidates_t;

GEM_INLINE void filtering_candidates_new(filtering_candidates_t* const filtering_candidates,mm_slab_t* const mm_slab);
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates);
GEM_INLINE void filtering_candidates_delete(filtering_candidates_t* const filtering_candidates);

GEM_INLINE void filtering_candidates_add_text_position(
    filtering_candidates_t* const filtering_candidates,const uint64_t index_position);
GEM_INLINE void filtering_candidates_add_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,const uint64_t interval_hi,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t region_errors);
GEM_INLINE void filtering_candidates_add_interval_set(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos);
GEM_INLINE void filtering_candidates_add_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t max_error);

GEM_INLINE uint64_t filtering_candidates_get_pending_candidates(filtering_candidates_t* const filtering_candidates);

/*
 * Filtering Candidates Verification
 */
GEM_INLINE void filtering_candidates_verify_pending(
    filtering_candidates_t* const filtering_candidates,matches_t* const matches,
    const bool use_levenshtein_distance,const bool store_hamming_matches);

#endif /* FILTERING_CANDIDATES_H_ */
