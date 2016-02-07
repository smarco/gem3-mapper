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
void gem_thread_register_id(const uint64_t thread_id) {
  MUTEX_BEGIN_SECTION(threads_id_locator_mutex) {
    if (threads_id_locator==NULL) threads_id_locator = ihash_new();
    const int64_t posix_thread_id = (int64_t)pthread_self();
    int64_t* gem_thread_id = ihash_get(threads_id_locator,posix_thread_id,int64_t);
    if (gem_thread_id==NULL) {
      gem_thread_id = mm_malloc(sizeof(int64_t));
      *gem_thread_id = thread_id;
      ihash_insert(threads_id_locator,posix_thread_id,gem_thread_id);
    } else {
      *gem_thread_id = thread_id;
    }
  } MUTEX_END_SECTION(threads_id_locator_mutex);
}
int64_t gem_thread_get_thread_id() {
  int64_t* const thread_id = ihash_get(threads_id_locator,(int64_t)pthread_self(),int64_t);
  return (thread_id) ? *thread_id : -1;
}
void gem_thread_cleanup() {
  if (threads_id_locator!=NULL) ihash_delete(threads_id_locator);
}
