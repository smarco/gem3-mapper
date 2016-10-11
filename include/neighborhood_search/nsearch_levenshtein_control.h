/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_LEVENSHTEIN_CONTROL_H_
#define NSEARCH_LEVENSHTEIN_CONTROL_H_

#include "utils/essentials.h"
#include "neighborhood_search/nsearch_schedule.h"
#include "neighborhood_search/nsearch_levenshtein_query.h"

/*
 * Search candidates cut-off
 */
bool nsearch_levenshtein_candidates_cutoff(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t num_candidates,
    nsearch_query_t* const nsearch_query,
    nsearch_query_t* const next_nsearch_query);

/*
 * Standard search terminate search-branch
 */
uint64_t nsearch_levenshtein_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t text_position,
    uint64_t lo,
    uint64_t hi,
    const uint64_t align_distance);

/*
 * Scheduled search terminate search-branch
 */
uint64_t nsearch_levenshtein_scheduled_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    nsearch_query_t* const nsearch_query,
    const uint64_t align_distance);

#endif /* NSEARCH_LEVENSHTEIN_CONTROL_H_ */
