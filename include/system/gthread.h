/*
 * PROJECT: GEMMapper
 * FILE: gthread.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
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
int64_t gem_thread_get_thread_id();
#define gtid() gem_thread_get_thread_id()

void gem_thread_cleanup();

/*
 * Mutex/ConditionalVars Helpers
 */
#define MUTEX_INIT(mutex) \
  gem_cond_fatal_error__perror(pthread_mutex_init(&(mutex),NULL),SYS_MUTEX_INIT);
#define MUTEX_DESTROY(mutex) \
  gem_cond_fatal_error__perror(pthread_mutex_destroy(&(mutex)),SYS_MUTEX_DESTROY);
#define MUTEX_BEGIN_SECTION(mutex) \
  gem_cond_fatal_error__perror(pthread_mutex_lock(&(mutex)),SYS_MUTEX);
#define MUTEX_END_SECTION(mutex) \
  gem_cond_fatal_error__perror(pthread_mutex_unlock(&(mutex)),SYS_MUTEX);

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
