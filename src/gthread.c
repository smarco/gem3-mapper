/*
 * PROJECT: GEMMapper
 * FILE: gthread.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "gthread.h"

#include "commons.h"
#include "mm.h"
#include "hash.h"

/*
 * Thread ID Locator
 */
ihash_t* threads_id_locator = NULL;
pthread_mutex_t threads_id_locator_mutex = PTHREAD_MUTEX_INITIALIZER;
inline void gem_thread_register_id(const uint64_t thread_id) {
  MUTEX_BEGIN_SECTION(threads_id_locator_mutex) {
    if (threads_id_locator==NULL) threads_id_locator = ihash_new();
    int64_t* const thread_id_int = mm_malloc_uint64();
    *thread_id_int = (int64_t)thread_id;
    ihash_insert(threads_id_locator,(int64_t)pthread_self(),thread_id_int);
  } MUTEX_END_SECTION(threads_id_locator_mutex);
}
inline int64_t gem_thread_get_thread_id() {
  int64_t* const thread_id = ihash_get(threads_id_locator,(int64_t)pthread_self(),int64_t);
  return (thread_id) ? *thread_id : -1;
}
