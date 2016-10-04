/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Sampled Run-Length (RL) data structure stores anchors every certain positions
 *   to enable translation between RL-space and text-space (original sequence before
 *   it was RL compacted)
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
