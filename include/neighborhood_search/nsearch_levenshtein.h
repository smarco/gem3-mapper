/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_LEVENSHTEIN_H_
#define NSEARCH_LEVENSHTEIN_H_

#include "utils/essentials.h"
#include "data_structures/interval_set.h"
#include "fm_index/fm_index.h"
#include "filtering/region_profile.h"
#include "neighborhood_search/dp_matrix.h"
#include "neighborhood_search/nsearch_schedule.h"

/*
 * Levenshtein Brute Force
 */
void nsearch_levenshtein_brute_force(
    approximate_search_t* const search,
    const bool supercondensed,
    matches_t* const matches);

/*
 * Levenshtein Neighborhood Search
 */
void nsearch_levenshtein(
    approximate_search_t* const search,
    matches_t* const matches);
void nsearch_levenshtein_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Display
 */
void nsearch_levenshtein_print_status(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches);

#endif /* NSEARCH_LEVENSHTEIN_H_ */
