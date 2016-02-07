/*
 * PROJECT: GEMMapper
 * FILE: region_profile_boost.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_PROFILE_BOOST_H_
#define REGION_PROFILE_BOOST_H_

#include "essentials.h"
#include "region_profile.h"

/*
 * Region Profile Generation
 */
void region_profile_generate_adaptive_boost(
    region_profile_t* const region_profile,fm_index_t* const fm_index,
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,mm_stack_t* const mm_stack);

#endif /* REGION_PROFILE_BOOST_H_ */
