/*
 * PROJECT: GEMMapper
 * FILE: sampled_sa.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef SAMPLED_SA_H_
#define SAMPLED_SA_H_

#include "essentials.h"

typedef enum { SAMPLING_RATE_1=0,  SAMPLING_RATE_2=1,   SAMPLING_RATE_4=2,
               SAMPLING_RATE_8=3,  SAMPLING_RATE_16=4,  SAMPLING_RATE_32=5,
               SAMPLING_RATE_64=6, SAMPLING_RATE_128=7, SAMPLING_RATE_256=8 } sampling_rate_t;

// Sampled-SA (Implementation independent)
typedef struct _sampled_sa_t sampled_sa_t;
typedef struct _sampled_sa_builder_t sampled_sa_builder_t;

/*
 * Loader
 */
GEM_INLINE sampled_sa_t* sampled_sa_read(fm_t* const file_manager);
GEM_INLINE sampled_sa_t* sampled_sa_read_mem(mm_t* const memory_manager);
GEM_INLINE void sampled_sa_delete(sampled_sa_t* const sampled_sa);

/*
 * Builder
 */
GEM_INLINE sampled_sa_builder_t* sampled_sa_builder_new(const uint64_t index_length,const sampling_rate_t sampling_rate);
GEM_INLINE uint64_t sampled_sa_builder_get_sampling_rate_value(sampled_sa_builder_t* const sampled_sa);
GEM_INLINE void sampled_sa_builder_set_sample(sampled_sa_builder_t* const sampled_sa,const uint64_t sa_sampled_position,const uint64_t sa_value);
GEM_INLINE void sampled_sa_builder_write(fm_t* const file_manager,sampled_sa_builder_t* const sampled_sa);
GEM_INLINE void sampled_sa_builder_delete(sampled_sa_builder_t* const sampled_sa);

/*
 * Accessors
 */
GEM_INLINE uint64_t sampled_sa_get_size(sampled_sa_t* const sampled_sa);
GEM_INLINE uint64_t sampled_sa_get_sampling_rate_value(sampled_sa_t* const sampled_sa);
GEM_INLINE sampling_rate_t sampled_sa_get_sampling_rate(sampled_sa_t* const sampled_sa);
GEM_INLINE bool sampled_sa_is_sampled(sampled_sa_t* const sampled_sa,const uint64_t sa_position);
GEM_INLINE uint64_t sampled_sa_get_sample(sampled_sa_t* const sampled_sa,const uint64_t sa_sampled_position);

/*
 * Display/Stats
 */
GEM_INLINE void sampled_sa_print(FILE* const stream,sampled_sa_t* const sampled_sa);
GEM_INLINE void sampled_sa_builder_print(FILE* const stream,sampled_sa_builder_t* const sampled_sa_builder);


#endif /* SAMPLED_SA_H_ */
