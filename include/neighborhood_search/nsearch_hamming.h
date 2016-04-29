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
#include "fm_index/fm_index.h"
#include "data_structures/interval_set.h"
#include "filtering/region_profile.h"
#include "neighborhood_search/nsearch_schedule.h"

/*
 * Hamming Brute Force
 */
void nsearch_hamming_brute_force(
    fm_index_t* const restrict fm_index,uint8_t* const restrict key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const restrict intervals_result,mm_stack_t* const restrict mm_stack);

/*
 * Perform scheduled search
 */
void nsearch_hamming_perform_scheduled_search(
    nsearch_schedule_t* const restrict nsearch_schedule,const uint64_t pending_searches,
    nsearch_operation_t* const restrict nsearch_operation,const uint64_t position,
    const uint64_t local_error,const uint64_t global_error);

/*
 * Hamming Neighborhood Search
 */
void nsearch_hamming(
    fm_index_t* const restrict fm_index,uint8_t* const restrict key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const restrict intervals_result,mm_stack_t* const restrict mm_stack);
void nsearch_hamming_preconditioned(
    fm_index_t* const restrict fm_index,region_profile_t* const restrict region_profile,
    uint8_t* const restrict key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const restrict intervals_result,mm_stack_t* const restrict mm_stack);



#endif /* NSEARCH_HAMMING_H_ */
