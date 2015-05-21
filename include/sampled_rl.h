/*
 * PROJECT: GEMMapper
 * FILE: sampled_rl.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef SAMPLED_RL_H_
#define SAMPLED_RL_H_

#include "essentials.h"

#include "sampled_sa.h"
#include "packed_integer_array.h"

/*
 * Checkers
 */
#define SAMPLED_RL_CHECK(sampled_rl) GEM_CHECK_NULL(sampled_rl)

/*
 * Constants
 */
#define SAMPLED_RL_MAX_RUN_LENGTH  5
#define SAMPLED_RL_SAMPLING_RATE  (1<<6)   /* 2^6=64 (HG => 460 MB)*/

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
GEM_INLINE sampled_rl_t* sampled_rl_new(
    const sampling_rate_t sampling_rate,const uint64_t num_samples,const uint64_t max_index);
GEM_INLINE sampled_rl_t* sampled_rl_read_mem(mm_t* const memory_manager);
GEM_INLINE void sampled_rl_write(fm_t* const file_manager,sampled_rl_t* const sampled_rl);
GEM_INLINE void sampled_rl_delete(sampled_rl_t* const sampled_rl);

/*
 * Accessors
 */
GEM_INLINE uint64_t sampled_rl_get_size(sampled_rl_t* const sampled_rl);
GEM_INLINE void sampled_rl_sample(sampled_rl_t* const sampled_rl,const uint64_t array_position,const uint64_t rl_position);
GEM_INLINE uint64_t sampled_rl_get_sample(sampled_rl_t* const sampled_rl,const uint64_t array_position);
/*
 * Display/Stats
 */
GEM_INLINE void sampled_rl_print(FILE* const stream,sampled_rl_t* const sampled_rl);

/*
 * Errors
 */
#define GEM_ERROR_SAMPLED_RL_WRONG_MODEL_NO "Sampled-RL error. Wrong Sampled-RL Model %lu (Expected model %lu)"

#endif /* SAMPLED_RL_H_ */
