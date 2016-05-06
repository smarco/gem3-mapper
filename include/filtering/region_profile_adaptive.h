/*
 * PROJECT: GEMMapper
 * FILE: region_profile_adaptive.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_PROFILE_ADAPTIVE_H_
#define REGION_PROFILE_ADAPTIVE_H_

#include "utils/essentials.h"
#include "filtering/region_profile.h"

/*
 * Region Profile Generator & Query
 */
typedef struct {
  // Region Profile
  region_profile_t* region_profile;
  // Region State
  uint64_t region_length;
  uint64_t last_cut;
  uint64_t lo_cut;
  uint64_t hi_cut;
  uint64_t expected_count;
  uint64_t max_steps;
  bool allow_zero_regions;
  // Query
  fm_index_t* fm_index;
  const uint8_t* key;
  uint64_t key_length;
  const bool* allowed_enc;
  // Query state
  uint64_t key_position;
  rank_mquery_t rank_mquery;
  uint64_t lo;
  uint64_t hi;
  // Mappabiliy
  double mappability_p_acc;
  uint64_t mappability_p_samples;
  double mappability_2p_acc;
  uint64_t mappability_2p_samples;
} region_profile_generator_t;

/*
 * Region Profile Adaptive Iterator
 */
void region_profile_generator_init(
    region_profile_generator_t* const generator,
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    const bool allow_zero_regions);
bool region_profile_generator_next_region(
    region_profile_t* const region_profile,
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model);

/*
 * Region Profile Generation
 */
void region_profile_generate_adaptive(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,
    const uint64_t max_regions,
    const bool allow_zero_regions);
void region_profile_generate_adaptive_limited(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    const region_profile_model_t* const profile_model,
    const uint64_t min_regions);

#endif /* REGION_PROFILE_ADAPTIVE_H_ */
