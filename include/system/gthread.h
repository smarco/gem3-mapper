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

#ifndef GTHREAD_H_
#define GTHREAD_H_

#include "system/commons.h"
#include <pthread.h>

/*
 * Types
 */
typedef void* (*pthread_handler_t)(void*);

/*
 * Thread ID Locator
 */
void gem_thread_register_id(const uint64_t thread_id);
int64_t gem_thread_get_thread_id(void);
#define gtid() gem_thread_get_thread_id()

void gem_thread_cleanup(void);

/*
 * Mutex/ConditionalVars Helpers
 */
#define MUTEX_INIT(mutex) \
  gem_cond_fatal_error(pthread_mutex_init(&(mutex),NULL),SYS_MUTEX_INIT);
#define MUTEX_DESTROY(mutex) \
  gem_cond_fatal_error(pthread_mutex_destroy(&(mutex)),SYS_MUTEX_DESTROY);
#define MUTEX_BEGIN_SECTION(mutex) \
  gem_cond_fatal_error(pthread_mutex_lock(&(mutex)),SYS_MUTEX);
#define MUTEX_END_SECTION(mutex) \
  gem_cond_fatal_error(pthread_mutex_unlock(&(mutex)),SYS_MUTEX);

#define CV_INIT(cv) \
  gem_cond_fatal_error(pthread_cond_init(&(cv),NULL),SYS_COND_VAR_INIT);
#define CV_DESTROY(cv) \
  gem_cond_fatal_error(pthread_cond_destroy(&(cv)),SYS_COND_VAR_INIT);
#define CV_SIGNAL(cv) \
  gem_cond_fatal_error(pthread_cond_signal(&(cv)),SYS_COND_VAR);
#define CV_BROADCAST(cv) \
  gem_cond_fatal_error(pthread_cond_broadcast(&(cv)),SYS_COND_VAR);
#define CV_WAIT(cv,mutex) \
  gem_cond_fatal_error(pthread_cond_wait(&(cv),&(mutex)),SYS_COND_VAR);

#endif /* GTHREAD_H_ */
