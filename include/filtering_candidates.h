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

#include "archive.h"
#include "locator.h"
#include "fm_index.h"
#include "pattern.h"
#include "bpm_align_gpu.h"
#include "approximate_search_parameters.h"
#include "region_profile.h"
#include "interval_set.h"
#include "filtering_region.h"
#include "matches.h"
#include "paired_matches.h"

/*
 * Filtering Candidates Vector
 */
typedef struct {
  /* Filtering Status */
  uint64_t total_candidates_accepted;         // Total number of candidates accepted
  /* Region Buffer */
  vector_t* regions_buffer;                   // Regions Buffer (region_t)
  /* Candidates */
  vector_t* filtering_positions;              // Candidate positions (filtering_position_t)
  ihash_t* verified_positions;                // Verified positions (uint64_t)
  vector_t* filtering_regions;                // Candidate regions (filtering_region_t)
  vector_t* verified_regions;                 // Verified regions (verified_region_t)
} filtering_candidates_t;

/*
 * Setup
 */
GEM_INLINE void filtering_candidates_init(filtering_candidates_t* const filtering_candidates);
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates);
GEM_INLINE void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates);

/*
 * Accessors
 */
GEM_INLINE uint64_t filtering_candidates_get_num_candidate_regions(const filtering_candidates_t* const filtering_candidates);
GEM_INLINE uint64_t filtering_candidates_count_candidate_regions(
    filtering_candidates_t* const filtering_candidates_end,const filtering_region_status_t filtering_region_status);

/*
 * Adding candidate positions
 */
GEM_INLINE void filtering_candidates_add_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,const uint64_t interval_hi,const uint64_t region_start_pos,
    const uint64_t region_end_pos,const uint64_t region_errors,mm_stack_t* const mm_stack);
GEM_INLINE void filtering_candidates_add_interval_set(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,mm_stack_t* const mm_stack);
GEM_INLINE void filtering_candidates_add_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t max_error,
    mm_stack_t* const mm_stack);

/*
 * Processing & Verification
 */
GEM_INLINE uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    mm_stack_t* const mm_stack);
GEM_INLINE uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack);

/*
 * BPM-Buffer API (Verification)
 */
GEM_INLINE uint64_t filtering_candidates_bpm_buffer_add(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,mm_stack_t* const mm_stack);
GEM_INLINE void filtering_candidates_bpm_buffer_align(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    pattern_t* const pattern,const strand_t search_strand,const search_actual_parameters_t* const search_actual_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end,
    matches_t* const matches,mm_stack_t* const mm_stack);

/*
 * Paired Verification
 */
GEM_INLINE void filtering_candidates_set_all_regions_pending(filtering_candidates_t* const filtering_candidates);
GEM_INLINE void filtering_candidates_set_all_regions_unverified(filtering_candidates_t* const filtering_candidates);
GEM_INLINE void filtering_candidates_paired_regions_filtering(
    filtering_candidates_t* const filtering_candidates_end1,filtering_candidates_t* const filtering_candidates_end2,
    const uint64_t min_template_length,const uint64_t max_template_length,const bool absolute_distance);
GEM_INLINE uint64_t filtering_candidates_extend_match(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    const match_trace_t* const extended_match,const pattern_t* const candidate_pattern,
    const strand_t candidate_search_strand,const bool search_onward,
    const search_actual_parameters_t* const candidate_actual_parameters,
    paired_matches_t* const paired_matches,const sequence_end_t candidate_end,
    mm_stack_t* const mm_stack);

/*
 * Display
 */
GEM_INLINE void filtering_candidates_print_matching_regions(
    FILE* const stream,filtering_candidates_t* const filtering_candidates);

#endif /* FILTERING_CANDIDATES_H_ */
