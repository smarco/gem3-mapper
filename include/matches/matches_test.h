/*
 * PROJECT: GEMMapper
 * FILE: matches.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_TEST_H_
#define MATCHES_TEST_H_

#include "utils/essentials.h"
#include "archive/search/archive_search_se_parameters.h"
#include "matches/matches.h"

/*
 * Matches Condition Tests
 */
bool matches_test_max_matches_reached(
    matches_t* const matches,
    const uint64_t mcs,
    const uint64_t key_length,
    search_parameters_t* const search_parameters);
bool matches_test_accuracy_reached(
    matches_t* const matches,
    const uint64_t mcs,
    const uint64_t key_length,
    search_parameters_t* const search_parameters,
    const uint64_t max_complete_error,
    uint64_t* const max_complete_error_required);

#endif /* MATCHES_TEST_H_ */
