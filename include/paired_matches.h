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
#include "matches.h"

/*
 * Checkers
 */
#define PAIRED_MATCHES_CHECK(paired_matches) GEM_CHECK_NULL(paired_matches)

/*
 * Paired Matches
 */
typedef struct {
  /* Search-matches state */
  uint64_t max_complete_stratum;            // MCS
  /* Text Collection Buffer */
  text_collection_t* text_collection;       // Stores text-traces (candidates/matches/regions/...)
  /* Matches Counters */
  vector_t* counters;                       // Global counters
  /* Single-End Matches */
  matches_t* matches_end1;                  // Matches end1
  matches_t* matches_end2;                  // Matches end2
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
GEM_INLINE matches_t* paired_matches_end1(paired_matches_t* const paired_matches);
GEM_INLINE matches_t* paired_matches_end2(paired_matches_t* const paired_matches);

/*
 * Counters
 */
GEM_INLINE void paired_matches_counters_add(
    paired_matches_t* const paired_matches,const uint64_t distance,const uint64_t num_matches);

/*
 * Adding Paired-Matches
 */
GEM_INLINE void paired_matches_add(
    paired_matches_t* const paired_matches,
    match_trace_t* const match_trace_end1,match_trace_t* const match_trace_end2);

/*
 * Error Messages
 */
//#define GEM_ERROR_PAIRED_MATCHES

#endif /* PAIRED_MATCHES_H_ */
