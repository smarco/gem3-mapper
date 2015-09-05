/*
 * PROJECT: GEMMapper
 * FILE: region_profile_fixed.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_PROFILE_FIXED_H_
#define REGION_PROFILE_FIXED_H_

#include "essentials.h"
#include "region_profile.h"

/*
 * Region Profile Schedule (generate the region partition)
 */
void region_profile_generate_fixed_schedule(
    region_profile_t* const region_profile,const uint8_t* const key,
    const uint64_t key_length,const bool* const allowed_enc,
    const uint64_t region_length);

/*
 * Region Profile Schedule (query the region partition into the index)
 */
void region_profile_generate_fixed_query(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key);

#endif /* REGION_PROFILE_FIXED_H_ */
