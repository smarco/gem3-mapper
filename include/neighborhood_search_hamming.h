/*
 * PROJECT: GEMMapper
 * FILE: neighborhood_search_hamming.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NEIGHBORHOOD_SEARCH_HAMMING_H_
#define NEIGHBORHOOD_SEARCH_HAMMING_H_

#include "essentials.h"
#include "fm_index.h"
#include "interval_set.h"

/*
 * Neighborhood Search
 */
void neighborhood_search_hamming_bidirectional(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack);

void neighborhood_search_hamming_brute_force(
    uint8_t* const key,const uint64_t key_length,const uint64_t max_error);

#endif /* NEIGHBORHOOD_SEARCH_HAMMING_H_ */
