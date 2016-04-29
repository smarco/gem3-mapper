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
#include "fm_index/fm_index.h"
#include "data_structures/interval_set.h"
#include "filtering/region_profile.h"
#include "neighborhood_search/nsearch_schedule.h"
#include "neighborhood_search/dp_matrix.h"

/*
 * Levenshtein Brute Force
 */
void nsearch_levenshtein_brute_force(
    fm_index_t* const restrict fm_index,uint8_t* const restrict key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const restrict intervals_result,mm_stack_t* const restrict mm_stack);

/*
 * Perform Levenshtein Scheduled Search
 */
uint64_t nsearch_levenshtein_perform_scheduled_search(
    nsearch_schedule_t* const restrict nsearch_schedule,const uint64_t pending_searches,
    nsearch_operation_t* const restrict nsearch_operation,const uint64_t global_error);

/*
 * Levenshtein Neighborhood Search
 */
void nsearch_levenshtein(
    fm_index_t* const restrict fm_index,uint8_t* const restrict key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const restrict intervals_result,mm_stack_t* const restrict mm_stack);
void nsearch_levenshtein_preconditioned(
    fm_index_t* const restrict fm_index,region_profile_t* const restrict region_profile,
    uint8_t* const restrict key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const restrict intervals_result,mm_stack_t* const restrict mm_stack);

/*
 * Display
 */
void nsearch_levenshtein_print_search_trace(
    FILE* const restrict stream,nsearch_schedule_t* const restrict nsearch_schedule,
    const uint64_t pending_searches);
void nsearch_levenshtein_print_pair_key_text(
    FILE* const restrict stream,nsearch_schedule_t* const restrict nsearch_schedule,
    nsearch_operation_t* const restrict nsearch_operation);


#endif /* NSEARCH_LEVENSHTEIN_H_ */
