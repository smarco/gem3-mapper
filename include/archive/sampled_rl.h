/*
 * PROJECT: GEMMapper
 * FILE: sampled_rl.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef SAMPLED_RL_H_
#define SAMPLED_RL_H_

#include "utils/essentials.h"
#include "utils/packed_integer_array.h"

/*
 * Constants
 */
#define SAMPLED_RL_SAMPLING_RATE   100

/*
 * Sampled Run-Length
 */
typedef struct {
  sampling_rate_t sampling_rate;                 // Sampling Rate
  packed_integer_array_t* packed_integer_array;  // Packed Samples
} sampled_rl_t;

/*
 * Loader/Setup
 */
sampled_rl_t* sampled_rl_new(
    const sampling_rate_t sampling_rate,
    const uint64_t num_samples,
    const uint64_t max_index);
sampled_rl_t* sampled_rl_read_mem(mm_t* const memory_manager);
void sampled_rl_write(
    fm_t* const file_manager,
    sampled_rl_t* const sampled_rl);
void sampled_rl_delete(sampled_rl_t* const sampled_rl);

/*
 * Accessors
 */
uint64_t sampled_rl_get_size(sampled_rl_t* const sampled_rl);
void sampled_rl_sample(
    sampled_rl_t* const sampled_rl,
    const uint64_t array_position,
    const uint64_t rl_position);
uint64_t sampled_rl_get_sample(
    sampled_rl_t* const sampled_rl,
    const uint64_t array_position);

/*
 * Display/Stats
 */
void sampled_rl_print(
    FILE* const stream,
    sampled_rl_t* const sampled_rl);

#endif /* SAMPLED_RL_H_ */
