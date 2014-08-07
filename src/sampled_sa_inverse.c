/*
 * PROJECT: GEMMapper
 * FILE: sampled_sa.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "sampled_sa.h"
#include "packed_integer_array.h"

/*
 * Checkers
 */
#define SAMPLED_SA_CHECK(sampled_sa) \
  GEM_CHECK_NULL(sampled_sa)

/*
 * Inverse-Sampled SA (Marked positions in SA space)
 */
struct _sampled_sa_t {
  uint64_t index_length;
  sampling_rate_t sampling_rate;
  packed_integer_array_t* packed_integer_array;
};
struct _sampled_sa_builder_t {
  // Meta-data
  uint64_t index_length;
  sampling_rate_t sampling_rate;
  // Integer Array (deflated)
  uint64_t num_sa_samples;
  uint64_t* raw_integer_array;
};

/*
 * Loader
 */
GEM_INLINE sampled_sa_t* sampled_sa_read(fm_t* const file_manager) {
  // Allocate handler
  sampled_sa_t* const sampled_sa = mm_alloc(sampled_sa_t);
  // Read Meta-Data
  sampled_sa->index_length = fm_read_uint64(file_manager);
  sampled_sa->sampling_rate = fm_read_uint64(file_manager);
  // Read PackedIntegerArray
  sampled_sa->packed_integer_array = packed_integer_array_read(file_manager);
  // Return
  return sampled_sa;
}
GEM_INLINE sampled_sa_t* sampled_sa_read_mem(mm_t* const memory_manager) {
  // Allocate handler
  sampled_sa_t* const sampled_sa = mm_alloc(sampled_sa_t);
  // Read Meta-Data
  sampled_sa->index_length = mm_read_uint64(memory_manager);
  sampled_sa->sampling_rate = mm_read_uint64(memory_manager);
  // Read PackedIntegerArray
  sampled_sa->packed_integer_array = packed_integer_array_read_mem(memory_manager);
  // Return
  return sampled_sa;
}
GEM_INLINE void sampled_sa_delete(sampled_sa_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  // Free PackedIntegerArray
  packed_integer_array_delete(sampled_sa->packed_integer_array);
  // Free handler
  mm_free(sampled_sa);
}
/*
 * Builder
 */
GEM_INLINE sampled_sa_builder_t* sampled_sa_builder_new(const uint64_t index_length,const sampling_rate_t sampling_rate) {
  // Allocate handler
  sampled_sa_builder_t* const sampled_sa = mm_alloc(sampled_sa_builder_t);
  sampled_sa->index_length = index_length;
  sampled_sa->sampling_rate = sampling_rate;
  // Allocate sampled SA
  sampled_sa->num_sa_samples = DIV_CEIL(index_length,1<<sampled_sa->sampling_rate);
  sampled_sa->raw_integer_array = mm_calloc(sampled_sa->num_sa_samples,uint64_t,false);
  // Return
  return sampled_sa;
}
GEM_INLINE uint64_t sampled_sa_builder_get_sampling_rate_value(sampled_sa_builder_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  return (1<<sampled_sa->sampling_rate);
}
GEM_INLINE void sampled_sa_builder_set_sample(sampled_sa_builder_t* const sampled_sa,const uint64_t sa_sampled_position,const uint64_t sa_value) {
  SAMPLED_SA_CHECK(sampled_sa);
  // Calculate array index position
  const uint64_t array_position = (sa_sampled_position >> sampled_sa->sampling_rate);
  // Store SA position
  sampled_sa->raw_integer_array[array_position] = sa_value;
}
GEM_INLINE void sampled_sa_builder_write(fm_t* const file_manager,sampled_sa_builder_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  // Write Meta-Data
  fm_write_uint64(file_manager,sampled_sa->index_length);
  fm_write_uint64(file_manager,sampled_sa->sampling_rate);
  // Write PackedIntegerArray
  packed_integer_array_write_adaptor(file_manager,
      sampled_sa->index_length,sampled_sa->num_sa_samples,sampled_sa->raw_integer_array);
}
GEM_INLINE void sampled_sa_builder_delete(sampled_sa_builder_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  mm_free(sampled_sa->raw_integer_array); // Free IntegerArray
  mm_free(sampled_sa); // Free handler
}
/*
 * Accessors
 */
GEM_INLINE uint64_t sampled_sa_get_size(const sampled_sa_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  return packed_integer_array_get_size(sampled_sa->packed_integer_array);
}
GEM_INLINE uint64_t sampled_sa_get_sampling_rate_value(const sampled_sa_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  return (1<<sampled_sa->sampling_rate);
}
GEM_INLINE sampling_rate_t sampled_sa_get_sampling_rate(const sampled_sa_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  return sampled_sa->sampling_rate;
}
GEM_INLINE bool sampled_sa_is_sampled(const sampled_sa_t* const sampled_sa,const uint64_t sa_position) {
  SAMPLED_SA_CHECK(sampled_sa);
  return ((sa_position % (1<<sampled_sa->sampling_rate)) == 0);
}
GEM_INLINE void sampled_sa_prefetch_sample(const sampled_sa_t* const sampled_sa,const uint64_t sa_sampled_position) {
  SAMPLED_SA_CHECK(sampled_sa);
  // Calculate array index position
  const uint64_t array_position = (sa_sampled_position>>sampled_sa->sampling_rate);
  // Prefetch SA position
  packed_integer_array_prefetch(sampled_sa->packed_integer_array,array_position);
}
GEM_INLINE uint64_t sampled_sa_get_sample(const sampled_sa_t* const sampled_sa,const uint64_t sa_sampled_position) {
  SAMPLED_SA_CHECK(sampled_sa);
  // Calculate array index position
  const uint64_t array_position = (sa_sampled_position>>sampled_sa->sampling_rate);
  // Store SA position
  return packed_integer_array_load(sampled_sa->packed_integer_array,array_position);
}
/*
 * Display/Stats
 */
GEM_INLINE void sampled_sa_print(FILE* const stream,sampled_sa_t* const sampled_sa) {
  SAMPLED_SA_CHECK(sampled_sa);
  tab_fprintf(stream,"[GEM]>Sampled-SA\n");
  tab_fprintf(stream,"  => Architecture sSA.pck\n");
  tab_fprintf(stream,"    => Sampled-Space SA \n");
  tab_fprintf(stream,"    => ArrayImpl     Packed-Array \n");
  tab_fprintf(stream,"  => SamplingRate 2^%lu\n",(uint64_t)sampled_sa->sampling_rate);
  tab_fprintf(stream,"  => Length %lu\n",sampled_sa->index_length);
  tab_fprintf(stream,"  => Size %lu MB\n",CONVERT_B_TO_MB(sampled_sa_get_size(sampled_sa)));
  tab_global_inc();
  packed_integer_array_print(stream,sampled_sa->packed_integer_array,false);
  tab_global_dec();
  // Flush
  fflush(stream);
}
GEM_INLINE void sampled_sa_builder_print(FILE* const stream,sampled_sa_builder_t* const sampled_sa_builder) {
  SAMPLED_SA_CHECK(sampled_sa_builder);
  tab_fprintf(stream,"[GEM]>Sampled-SA\n");
  tab_fprintf(stream,"  => Architecture sSA.raw\n");
  tab_fprintf(stream,"    => Sampled-Space SA \n");
  tab_fprintf(stream,"    => ArrayImpl     Raw-Array \n");
  tab_fprintf(stream,"  => SamplingRate 2^%lu\n",(uint64_t)sampled_sa_builder->sampling_rate);
  tab_fprintf(stream,"  => Length %lu\n",sampled_sa_builder->index_length);
  tab_fprintf(stream,"  => Size %lu MB\n",CONVERT_B_TO_MB(sampled_sa_builder->index_length*UINT64_SIZE));
  // Flush
  fflush(stream);
}

