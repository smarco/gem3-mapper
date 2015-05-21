/*
 * PROJECT: GEMMapper
 * FILE: paired_matches.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PAIRED_MATCHES_H_
#define PAIRED_MATCHES_H_

#include "essentials.h"
#include "interval_set.h"
#include "text_collection.h"
#include "mapper_stats.h"
#include "matches.h"
#include "approximate_search_parameters.h"

/*
 * Checkers
 */
#define PAIRED_MATCHES_CHECK(paired_matches) GEM_CHECK_NULL(paired_matches)

/*
 * Paired Matches
 */
typedef struct {
  match_trace_t* match_end1;           // Map end1
  match_trace_t* match_end2;           // Map end2
  uint8_t mapq_score;                  // MAPQ score
  uint64_t template_length;            // Template observed length
  uint64_t distance;                   // Distance of the paired-alignment
  pair_orientation_t pair_orientation; // Pair orientation (concordant/discordant)
} paired_match_t;
typedef struct {
  /* State */
  uint64_t max_complete_stratum;            // MCS
  /* Text Collection Buffer */
  text_collection_t* text_collection;       // Stores text-traces (candidates/matches/regions/...)
  /* Matches Counters */
  vector_t* counters;                       // Global counters
  vector_t* discordant_counters;            // Auxiliary counters for discordant matches
  /* Single-End Matches */
  matches_t* matches_end1;                  // Matches end1
  matches_t* matches_end2;                  // Matches end2
  /* Paired-End Matches */
  vector_t* matches;                        // (paired_match_t)
  vector_t* discordant_matches;             // Auxiliary for discordant matches (paired_match_t)
} paired_matches_t;

/*
 * Setup
 */
GEM_INLINE paired_matches_t* paired_matches_new();
GEM_INLINE void paired_matches_configure(paired_matches_t* const paired_matches,text_collection_t* const text_collection);
GEM_INLINE void paired_matches_clear(paired_matches_t* const paired_matches);
GEM_INLINE void paired_matches_delete(paired_matches_t* const paired_matches);

/*
 * Accessors
 */
GEM_INLINE bool paired_matches_is_mapped(paired_matches_t* const paired_matches);

/*
 * Adding Paired-Matches
 */
GEM_INLINE void paired_matches_add(
    paired_matches_t* const paired_matches,match_trace_t* const match_trace_end1,
    match_trace_t* const match_trace_end2,const pair_orientation_t pair_orientation,const uint64_t template_length);

/*
 * Finding Pairs
 */
GEM_INLINE void paired_matches_pair_match_with_mates(
    paired_matches_t* const paired_matches,search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,const pair_orientation_t pair_orientation,
    match_trace_t* const match_trace,const sequence_end_t mate_end,
    match_trace_t* const mates_array,const uint64_t num_mates_trace);
GEM_INLINE void paired_matches_find_pairs(
    paired_matches_t* const paired_matches,search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats);
GEM_INLINE void paired_matches_find_discordant_pairs(
    paired_matches_t* const paired_matches,search_parameters_t* const search_parameters);

/*
 * Utils
 */
GEM_INLINE void paired_matches_sort_by_distance(paired_matches_t* const paired_matches);
GEM_INLINE void paired_matches_sort_by_mapq_score(paired_matches_t* const paired_matches);

GEM_INLINE uint64_t paired_match_get_template_observed_length(
    const uint64_t begin_position_1,const uint64_t end_position_1,
    const uint64_t begin_position_2,const uint64_t end_position_2);

GEM_INLINE uint64_t paired_match_calculate_distance(
    const match_trace_t* const match_trace_end1,const match_trace_t* const match_trace_end2);

/*
 * Error Messages
 */
//#define GEM_ERROR_PAIRED_MATCHES

#endif /* PAIRED_MATCHES_H_ */
