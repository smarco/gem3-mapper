/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "nsearch_hamming.h"
#include "nsearch_partition.h"
#include "nsearch_schedule.h"

/*
 * Brute Force
 */
GEM_INLINE void nsearch_levenshtein_brute_force_search(
    char* const search_string,uint8_t* const key,
    const uint64_t key_length,const uint64_t current_position,
    const uint64_t current_error,const uint64_t max_error) {
}
GEM_INLINE void nsearch_levenshtein_brute_force(uint8_t* const key,const uint64_t key_length,const uint64_t max_error) {
}
/*
 * Perform scheduled search
 */
GEM_INLINE void nsearch_levenshtein_enumerate_perform_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule,const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,const uint64_t position,
    const uint64_t local_error,const uint64_t global_error) {
}
/*
 * Hamming Neighborhood Search
 */
GEM_INLINE void nsearch_levenshtein(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
}
/*
 * Neighborhood Search (Preconditioned by region profile)
 */
GEM_INLINE void nsearch_levenshtein_preconditioned(
    fm_index_t* const fm_index,region_profile_t* const region_profile,
    uint8_t* const key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
}
