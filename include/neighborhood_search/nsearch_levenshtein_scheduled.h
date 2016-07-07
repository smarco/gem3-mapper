/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein_scheduled.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_LEVENSHTEIN_SCHEDULED_H_
#define NSEARCH_LEVENSHTEIN_SCHEDULED_H_

#include "utils/essentials.h"
#include "data_structures/interval_set.h"
#include "fm_index/fm_index.h"
#include "fm_index/fm_index_query.h"
#include "filtering/region_profile.h"
#include "neighborhood_search/dp_matrix.h"
#include "neighborhood_search/nsearch_schedule.h"

/*
 * Levenshtein Scheduled Search Query
 */
uint64_t nsearch_levenshtein_scheduled_search_operation_query(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    fm_2interval_t* const fm_2interval);

/*
 * Levenshtein Scheduled Search
 */
uint64_t nsearch_levenshtein_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const current_nsearch_operation,
    fm_2interval_t* const fm_2interval);

#endif /* NSEARCH_LEVENSHTEIN_SCHEDULED_H_ */
