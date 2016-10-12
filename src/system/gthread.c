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
 */

#include "system/gthread.h"
#include "system/mm.h"
#include "utils/hash.h"

/*
 * Thread ID Locator
 */
ihash_t* threads_id_locator = NULL;
pthread_mutex_t threads_id_locator_mutex = PTHREAD_MUTEX_INITIALIZER;
void gem_thread_register_id(const uint64_t thread_id) {
  MUTEX_BEGIN_SECTION(threads_id_locator_mutex) {
    if (threads_id_locator==NULL) threads_id_locator = ihash_new(NULL);
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
int64_t gem_thread_get_thread_id(void) {
  int64_t* const thread_id = ihash_get(threads_id_locator,(int64_t)pthread_self(),int64_t);
  return (thread_id) ? *thread_id : -1;
}
void gem_thread_cleanup(void) {
  if (threads_id_locator!=NULL) {
    IHASH_BEGIN_ITERATE(threads_id_locator,it_key,it_element,int64_t) {
      mm_free(it_element);
    } IHASH_END_ITERATE;
    ihash_delete(threads_id_locator);
  }
}
