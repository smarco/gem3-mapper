/*
 * PROJECT: GEMMapper
 * FILE: nsearch_hamming.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_HAMMING_H_
#define NSEARCH_HAMMING_H_

#include "utils/essentials.h"
#include "data_structures/interval_set.h"
#include "fm_index/fm_index.h"
#include "filtering/region_profile.h"
#include "neighborhood_search/nsearch_schedule.h"

/*
 * Hamming Brute Force
 */
void nsearch_hamming_brute_force(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Perform scheduled search
 */
void nsearch_hamming_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule);

/*
 * Hamming Neighborhood Search
 */
void nsearch_hamming(
    approximate_search_t* const search,
    matches_t* const matches);
void nsearch_hamming_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* NSEARCH_HAMMING_H_ */
