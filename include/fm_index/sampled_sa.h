/*
 * PROJECT: GEMMapper
 * FILE: sampled_sa.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
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
    sampled_sa_t* const sampled_sa,
    const bool display_data);
void sampled_sa_builder_print(
    FILE* const stream,
    sampled_sa_builder_t* const sampled_sa);

#endif /* SAMPLED_SA_H_ */
