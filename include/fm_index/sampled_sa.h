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
 *   Implements a data structure that stores SA samples each certain
 *   sampling-rate. Samples are stored bit-compacted
 */

#ifndef SAMPLED_SA_H_
#define SAMPLED_SA_H_

#include "utils/essentials.h"
#include "utils/packed_integer_array.h"

/*
 * Sampled-SA
 */
typedef struct {
  // Meta-data
  uint64_t index_length;
  sampling_rate_t sa_sampling_rate;
  sampling_rate_t text_sampling_rate;
  // Integer Array
  packed_integer_array_t* packed_integer_array;
} sampled_sa_t;

/*
 * Sampled-SA Builder
 */
typedef struct {
  // Meta-data
  uint64_t index_length;
  sampling_rate_t sa_sampling_rate;
  sampling_rate_t text_sampling_rate;
  // Integer Array Builder
  uint64_t num_chunks;
  packed_integer_array_builder_t** array_builder;
  // Sampled Positions Bitmap
  pthread_mutex_t sampled_bitmap_mutex;
  uint64_t sampled_bitmap_size;
  uint64_t* sampled_bitmap_mem;
  // Raw Samples (SA-Space Samples)
  uint64_t* sa_raw_samples;
  uint64_t sa_raw_samples_length;
} sampled_sa_builder_t;

/*
 * Loader
 */
sampled_sa_t* sampled_sa_read_mem(mm_t* const memory_manager);
void sampled_sa_write(fm_t* const file_manager,sampled_sa_t* const sampled_sa);
void sampled_sa_delete(sampled_sa_t* const sampled_sa);

/*
 * Builder
 */
sampled_sa_builder_t* sampled_sa_builder_new(
    const uint64_t index_length,
    const uint64_t num_builders,
    const sampling_rate_t sa_sampling_rate,
    const sampling_rate_t text_sampling_rate,
    const bool generate_raw_sampled_sa,
    mm_slab_t* const mm_slab);
void sampled_sa_builder_delete_samples(sampled_sa_builder_t* const sampled_sa);
void sampled_sa_builder_delete(sampled_sa_builder_t* const sampled_sa);

void sampled_sa_builder_set_sample(
    sampled_sa_builder_t* const sampled_sa,
    const uint64_t chunk_number,
    const uint64_t array_position,
    const uint64_t sa_value);
void sampled_sa_builder_write(
    fm_t* const file_manager,
    sampled_sa_builder_t* const sampled_sa);

uint64_t sampled_sa_builder_get_sa_sampling_rate(const sampled_sa_builder_t* const sampled_sa);
uint64_t sampled_sa_builder_get_text_sampling_rate(const sampled_sa_builder_t* const sampled_sa);
uint64_t* sampled_sa_builder_get_sampled_bitmap(const sampled_sa_builder_t* const sampled_sa);

/*
 * Accessors
 */
uint64_t sampled_sa_get_size(const sampled_sa_t* const sampled_sa);
uint64_t sampled_sa_get_sa_sampling_rate(const sampled_sa_t* const sampled_sa);
uint64_t sampled_sa_get_text_sampling_rate(const sampled_sa_t* const sampled_sa);

void sampled_sa_prefetch_sample(
    const sampled_sa_t* const sampled_sa,
    const uint64_t array_position);
uint64_t sampled_sa_get_sample(
    const sampled_sa_t* const sampled_sa,
    const uint64_t array_position);
void sampled_sa_set_sample(
    sampled_sa_t* const sampled_sa,
    const uint64_t array_position,
    const uint64_t sa_value);

/*
 * Display/Stats
 */
void sampled_sa_print(
    FILE* const stream,
    sampled_sa_t* const sampled_sa);
void sampled_sa_builder_print(
    FILE* const stream,
    sampled_sa_builder_t* const sampled_sa);

#endif /* SAMPLED_SA_H_ */
