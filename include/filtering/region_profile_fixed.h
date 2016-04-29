/*
 * PROJECT: GEMMapper
 * FILE: region_profile_fixed.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_PROFILE_FIXED_H_
#define REGION_PROFILE_FIXED_H_

#include "utils/essentials.h"
#include "filtering/region_profile.h"

/*
 * Region Profile Schedule (generate the region partition)
 */
void region_profile_generate_fixed_partition(
    region_profile_t* const restrict region_profile,
    const uint8_t* const restrict key,
    const uint64_t key_length,
    const bool* const restrict allowed_enc,
    const uint64_t min_region_length);

/*
 * Region Profile Schedule (query the region partition into the index)
 */
void region_profile_generate_fixed_query(
    region_profile_t* const restrict region_profile,
    fm_index_t* const restrict fm_index,
    const uint8_t* const restrict key);

/*
 * Display/Benchmark
 */
void region_profile_print_mappability(
    FILE* const restrict stream,
    fm_index_t* const restrict fm_index,
    const bool* const restrict allowed_enc,
    const uint8_t* key,
    const uint64_t key_length,
    const bool print_profiles,
    mm_stack_t* const restrict mm_stack);
void region_profile_print_benchmark(
    FILE* const restrict stream,
    const region_profile_t* const restrict region_profile,
    fm_index_t* const restrict fm_index,
    const uint8_t* key);

#endif /* REGION_PROFILE_FIXED_H_ */
