/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_LEVENSHTEIN_QUERY_H_
#define NSEARCH_LEVENSHTEIN_QUERY_H_

#include "utils/essentials.h"
#include "fm_index/fm_index_query.h"
#include "neighborhood_search/nsearch_schedule.h"

/*
 * NSearch query helper
 */
typedef struct {
  fm_2interval_t fm_2interval;
  fm_2erank_elms_t lo_2erank_elms;
  fm_2erank_elms_t hi_2erank_elms;
  uint64_t max_steps;
} nsearch_query_t;

/*
 * Standard search query
 */
void nsearch_levenshtein_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t current_position,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out);

/*
 * Scheduled search precomputed-query
 */
void nsearch_levenshtein_scheduled_precompute_query(
    nsearch_schedule_t* const nsearch_schedule,
    const bool forward_search,
    nsearch_query_t* const nsearch_query);

/*
 * Scheduled search query
 */
uint64_t nsearch_levenshtein_scheduled_query_exact(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    const uint8_t char_enc,
    nsearch_query_t* const nsearch_query);
uint64_t nsearch_levenshtein_scheduled_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    const uint8_t char_enc,
    nsearch_query_t* const nsearch_query,
    nsearch_query_t* const next_nsearch_query);

#endif /* NSEARCH_LEVENSHTEIN_QUERY_H_ */
