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

#include "system/mm_pool.h"
#include "system/commons.h"
#include "system/errors.h"
#include "system/fm.h"
#include "system/gthread.h"
#include "utils/vector.h"

/*
 * MM-Pool
 */
typedef struct {
  /* Huge Slab (128MB pages in 512MB)*/
  mm_slab_t* mm_slab_128MB;
  /* Extra-Large Slab (32MB pages in 256MB) */
  mm_slab_t* mm_slab_32MB;
  /* Large Slab (8MB pages in 64MB) */
  mm_slab_t* mm_slab_8MB;
  /* Regular Slab (2MB pages in 32MB) */
  mm_slab_t* mm_slab_2MB;
  /* Small Slab (128KB pages in 32MB)*/
  mm_slab_t* mm_slab_128KB;
  /* Mutex */
  pthread_mutex_t mutex;
} mm_pool_t;

/*
 * Global memory pool
 */
mm_pool_t* gem_memory_pool = NULL;

/*
 * Setup
 */
mm_slab_t* mm_pool_get_slab(const mm_pool_type_t mm_pool_type) {
  // Check @gem_memory_pool
  if (gem_memory_pool==NULL) {
    gem_memory_pool = mm_alloc(mm_pool_t);
    gem_memory_pool->mm_slab_128MB = NULL;
    gem_memory_pool->mm_slab_32MB = NULL;
    gem_memory_pool->mm_slab_8MB = NULL;
    gem_memory_pool->mm_slab_2MB = NULL;
    gem_memory_pool->mm_slab_128KB = NULL;
    MUTEX_INIT(gem_memory_pool->mutex);
  }
  // Critical Section
  mm_slab_t* mm_slab;
  MUTEX_BEGIN_SECTION(gem_memory_pool->mutex) {
    // Choose slab
    switch (mm_pool_type) {
      case mm_pool_128KB: // Small
        if (gem_memory_pool->mm_slab_128KB==NULL) {
          gem_memory_pool->mm_slab_128KB = mm_slab_new_(BUFFER_SIZE_128K,BUFFER_SIZE_32M,MM_UNLIMITED_MEM,"small.128KB");
        }
        mm_slab = gem_memory_pool->mm_slab_128KB;
        break;
      case mm_pool_2MB: // Regular
        if (gem_memory_pool->mm_slab_2MB==NULL) {
          gem_memory_pool->mm_slab_2MB = mm_slab_new_(BUFFER_SIZE_2M,BUFFER_SIZE_32M,MM_UNLIMITED_MEM,"regular.2MB");
        }
        mm_slab = gem_memory_pool->mm_slab_2MB;
        break;
      case mm_pool_8MB: // Large
        if (gem_memory_pool->mm_slab_8MB==NULL) {
          gem_memory_pool->mm_slab_8MB = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_64M,MM_UNLIMITED_MEM,"large.8MB");
        }
        mm_slab = gem_memory_pool->mm_slab_8MB;
        break;
      case mm_pool_32MB: // Extra-large
        if (gem_memory_pool->mm_slab_32MB==NULL) {
          gem_memory_pool->mm_slab_32MB = mm_slab_new_(BUFFER_SIZE_32M,BUFFER_SIZE_256M,MM_UNLIMITED_MEM,"huge.32MB");
        }
        mm_slab = gem_memory_pool->mm_slab_32MB;
        break;
      case mm_pool_128MB: // Huge
        if (gem_memory_pool->mm_slab_128MB==NULL) {
          gem_memory_pool->mm_slab_128MB = mm_slab_new_(BUFFER_SIZE_128M,BUFFER_SIZE_256M,MM_UNLIMITED_MEM,"huge.128MB");
        }
        mm_slab = gem_memory_pool->mm_slab_128MB;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  } MUTEX_END_SECTION(gem_memory_pool->mutex);
  return mm_slab;
}
void mm_pool_delete() {
  // Check @gem_memory_pool
  if (gem_memory_pool!=NULL) {
    // Delete slabs
    if (gem_memory_pool->mm_slab_128KB!=NULL) mm_slab_delete(gem_memory_pool->mm_slab_128KB);
    if (gem_memory_pool->mm_slab_2MB!=NULL) mm_slab_delete(gem_memory_pool->mm_slab_2MB);
    if (gem_memory_pool->mm_slab_8MB!=NULL) mm_slab_delete(gem_memory_pool->mm_slab_8MB);
    if (gem_memory_pool->mm_slab_32MB!=NULL) mm_slab_delete(gem_memory_pool->mm_slab_32MB);
    if (gem_memory_pool->mm_slab_128MB!=NULL) mm_slab_delete(gem_memory_pool->mm_slab_128MB);
    // Destroy mutex
    MUTEX_DESTROY(gem_memory_pool->mutex);
    // Free handler
    mm_free(gem_memory_pool);
  }
}

