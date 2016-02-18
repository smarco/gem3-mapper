/*
 * PROJECT: GEMMapper
 * FILE: mm_pool.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *     - PoolMemory
 *         Pool of Slabs as to gather all slabs needed along a program (of different segment size)
 *         This structure (essentially a hub) concentrates the memory sources of a program trying to avoid
 *           the individual usage of slabs and increase the fragmentation and unused blocks.
 *         Many slab sizes are available { mm_pool_128KB , mm_pool_2MB , mm_pool_8MB, mm_pool_32MB, mm_pool_128MB }.
 *           If bigger blocks are required, tailored mm_slabs should be created for that purpose
 *           (as no possible reuse of those slab can be expected)
 */

#ifndef MM_POOL_H_
#define MM_POOL_H_

#include "system/mm_slab.h"

typedef enum { mm_pool_128KB , mm_pool_2MB , mm_pool_8MB, mm_pool_32MB, mm_pool_128MB } mm_pool_type_t;

/*
 * Setup
 */
mm_slab_t* mm_pool_get_slab(const mm_pool_type_t mm_pool_type);
void mm_pool_delete();

#endif /* MM_POOL_H_ */
