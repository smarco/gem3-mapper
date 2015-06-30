/*
 * PROJECT: GEMMapper
 * FILE: paired_map.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PAIRED_MAP_H_
#define PAIRED_MAP_H_

#include "essentials.h"
#include "matches.h"
#include "paired_search_parameters.h"

/*
 * Paired Matches
 */
typedef struct {
  uint64_t match_end1_offset;          // Map end1
  uint64_t match_end2_offset;          // Map end2
  uint8_t mapq_score;                  // MAPQ score
  uint64_t template_length;            // Template observed length
  double template_length_sigma;        // Number of sigmas deviation from the std-distribution
  uint64_t distance;                   // Distance of the paired-alignment
  uint64_t edit_distance;              // Edit Distance of the paired-alignment
  int32_t swg_score;                   // Distance of the paired-alignment
  pair_orientation_t pair_orientation; // Pair orientation (concordant/discordant)
} paired_map_t;

/*
 * Accessors
 */
GEM_INLINE uint64_t paired_map_compute_distance(
    const match_trace_t* const match_end1,const match_trace_t* const match_end2);
GEM_INLINE uint64_t paired_map_compute_edit_distance(
    const match_trace_t* const match_end1,const match_trace_t* const match_end2);
GEM_INLINE uint64_t paired_map_compute_swg_score(
    const match_trace_t* const match_end1,const match_trace_t* const match_end2);

GEM_INLINE uint64_t paired_map_get_distance(paired_map_t* const paired_map);
GEM_INLINE uint64_t paired_map_get_edit_distance(paired_map_t* const paired_map);
GEM_INLINE int32_t paired_map_get_swg_score(paired_map_t* const paired_map);

/*
 * Utils
 */
GEM_INLINE uint64_t paired_map_get_template_observed_length(
    const uint64_t begin_position_1,const uint64_t end_position_1,
    const uint64_t begin_position_2,const uint64_t end_position_2);

#endif /* PAIRED_MAP_H_ */
