/*
 * PROJECT: GEMMapper
 * FILE: paired_map.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PAIRED_MAP_H_
#define PAIRED_MAP_H_

#include "utils/essentials.h"
#include "archive/archive_search_paired_parameters.h"
#include "matches/matches.h"

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
  pair_relation_t pair_relation;       // Pair relation (concordant/discordant)
  pair_orientation_t pair_orientation; // Pair orientation (FR,RF,FF,RR)
  pair_layout_t pair_layout;           // Pair layout (pair_layout_separate,pair_layout_overlap,pair_layout_contain)
  uint64_t index_position;             // Pair index position (only used for sorting purposes)
} paired_map_t;

/*
 * Accessors
 */
uint64_t paired_map_compute_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2);
uint64_t paired_map_compute_edit_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2);
uint64_t paired_map_compute_swg_score(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2);

uint64_t paired_map_get_distance(paired_map_t* const paired_map);
uint64_t paired_map_get_edit_distance(paired_map_t* const paired_map);
int32_t paired_map_get_swg_score(paired_map_t* const paired_map);

#endif /* PAIRED_MAP_H_ */
