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

#include "locator.h"
#include "fm_index.h"
#include "dna_text.h"
#include "pattern.h"
#include "bpm_align_gpu.h"

#include "approximate_search_parameters.h"
#include "region_profile.h"
#include "interval_set.h"
#include "matches.h"

/*
 * Filtering Candidates Vector
 */
typedef struct {
  /* Pending candidates */
  vector_t* candidate_positions;              // Candidates positions (candidate_position_t)
  /* Checked Positions */
  vector_t* verified_candidate_positions;     // Verified positions (uint64_t)
  /* Candidates accepted */
  uint64_t num_candidates_accepted;           // Total number of candidates accepted
  uint64_t max_candidates_accepted;           // Quick abandon due to maximum candidates accepted
  /* Internals */
  vector_t* regions_buffer;                   // Regions Buffer (region_t)
  text_collection_t* candidates_collection;   // Candidates Text-Collection (Stores candidates Texts)
} filtering_candidates_t;

GEM_INLINE void filtering_candidates_new(filtering_candidates_t* const filtering_candidates);
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates);
GEM_INLINE void filtering_candidates_delete(filtering_candidates_t* const filtering_candidates);

GEM_INLINE void filtering_candidates_set_max_candidates_accepted(
    filtering_candidates_t* const filtering_candidates,const uint64_t max_candidates_accepted);
GEM_INLINE bool filtering_candidates_is_max_candidates_reached(filtering_candidates_t* const filtering_candidates);

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
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,
    const approximate_search_parameters_t* const search_parameters,matches_t* const matches);
GEM_INLINE uint64_t filtering_candidates_add_to_bpm_buffer(
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,
    const approximate_search_parameters_t* const search_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer);

#endif /* FILTERING_CANDIDATES_H_ */
