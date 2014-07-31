/*
 * PROJECT: GEMMapper
 * FILE: profiler.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple time/functional profiler module
 */

#include "profiler.h"
#include "gthread.h"

/*
 * Timers
 */
GEM_INLINE void START_TIMER(gem_timer_t* const timer) {
  gettimeofday(&(timer->begin_timer),NULL);
}
GEM_INLINE void STOP_TIMER(gem_timer_t* const timer) {
  gettimeofday(&(timer->end_timer),NULL);
  timer->time_acc += TIME_DIFF(timer->begin_timer,timer->end_timer);
}
GEM_INLINE void RESET_TIMER(gem_timer_t* const timer) {
  timer->time_acc = 0;
}
GEM_INLINE void RESET__START_TIMER(gem_timer_t* const timer) {
  RESET_TIMER(timer);
  START_TIMER(timer);
}
GEM_INLINE double GET_TIMER(gem_timer_t* const timer) {
  return timer->time_acc;
}

#ifndef GEM_NOPROFILE
/*
 * Constants
 */
#define GP_MAX_COUNTERS 1000

/*
 * Profile
 */
typedef struct {
  // Time counters
  struct timeval* begin_timers;
  struct timeval* end_timers;
  double* timers;
  // General counters & functions
  uint64_t* counters;
} profile_t;
typedef struct {
  // Profiler
  profile_t* profile;
  // Limits
  uint64_t num_threads;
} profiler_t;

// THE GREAT PROFILER
profiler_t gem_profile;

/*
 * Setup
 */
GEM_INLINE void PROF_NEW(const uint64_t num_threads) {
  // Allocate handler
  gem_profile.profile = mm_calloc(num_threads,profile_t,true);
  // Initialize
  gem_profile.num_threads = num_threads;
  // Allocate all profiles
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    gem_profile.profile[i].begin_timers = mm_calloc(GP_MAX_COUNTERS,struct timeval,true);
    gem_profile.profile[i].end_timers = mm_calloc(GP_MAX_COUNTERS,struct timeval,true);
    gem_profile.profile[i].timers = mm_calloc(GP_MAX_COUNTERS,double,true);
    gem_profile.profile[i].counters = mm_calloc(GP_MAX_COUNTERS,uint64_t,true);
  }
}
GEM_INLINE void PROF_DELETE() {
  // Release all profile counters
  uint64_t i;
  for (i=0;i<gem_profile.num_threads;++i) {
    mm_free(gem_profile.profile[i].begin_timers);
    mm_free(gem_profile.profile[i].end_timers);
    mm_free(gem_profile.profile[i].timers);
    mm_free(gem_profile.profile[i].counters);
  }
  mm_free(gem_profile.profile);
}
/*
 * TIME functions
 */
GEM_INLINE void PROF_START_TIMER(const uint64_t timer) {
  gettimeofday((gem_profile.profile[gtid()].begin_timers+timer),NULL);
}
GEM_INLINE void PROF_STOP_TIMER(const uint64_t timer) {
  const uint64_t gtid = gtid();
  gettimeofday((gem_profile.profile[gtid].end_timers+timer),NULL);
  gem_profile.profile[gtid].timers[timer] +=
      TIME_DIFF(gem_profile.profile[gtid].begin_timers[timer],gem_profile.profile[gtid].end_timers[timer]);
}
GEM_INLINE void PROF_RESET_TIMER(const uint64_t timer) {
  gem_profile.profile[gtid()].timers[timer]=0.0;
}
GEM_INLINE double PROF_GET_TIMER(const uint64_t timer) {
  return gem_profile.profile[gtid()].timers[timer];
}
/*
 * GENERAL COUNTERS functions
 */
GEM_INLINE void PROF_RESET_COUNTER(const uint64_t counter) {
  gem_profile.profile[gtid()].counters[counter] = 0;
}
GEM_INLINE uint64_t PROF_GET_COUNTER(const uint64_t counter) {
  return gem_profile.profile[gtid()].counters[counter];
}
GEM_INLINE void PROF_SET_COUNTER(const uint64_t counter,const uint64_t value) {
  gem_profile.profile[gtid()].counters[counter] = value;
}
GEM_INLINE void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value) {
  gem_profile.profile[gtid()].counters[counter] += value;
}
GEM_INLINE void PROF_INC_COUNTER(const uint64_t counter) {
  ++(gem_profile.profile[gtid()].counters[counter]);
}
GEM_INLINE void PROF_DEC_COUNTER(const uint64_t counter) {
  --(gem_profile.profile[gtid()].counters[counter]);
}
/*
 * AGGREGATED functions
 */
GEM_INLINE void PROF_BEGIN(const uint64_t prof_counter) {
  PROF_INC_COUNTER(prof_counter);
  PROF_START_TIMER(prof_counter);
}
GEM_INLINE void PROF_END(const uint64_t prof_counter) {
  PROF_STOP_TIMER(prof_counter);
}
/*
 * Display statistics
 */
GEM_INLINE double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t tota_timer) {
  const uint64_t gtid = gtid();
  return PERCENTAGE(gem_profile.profile[gtid].timers[timer],gem_profile.profile[gtid].timers[tota_timer]);
}
GEM_INLINE uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter) {
  const uint64_t gtid = gtid();
  return PERCENTAGE(gem_profile.profile[gtid].counters[counter],gem_profile.profile[gtid].counters[total_counter]);
}
GEM_INLINE double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2) {
  const uint64_t gtid = gtid();
  return (float)gem_profile.profile[gtid].counters[counter1]/(float)gem_profile.profile[gtid].counters[counter2];
}
GEM_INLINE double PROF_TIME_PER_CALL(const uint64_t timer,const uint64_t counter) {
  const uint64_t gtid = gtid();
  return (gem_profile.profile[gtid].counters[counter]!=0) ?
      1000.0*((double)gem_profile.profile[gtid].timers[timer] /
              (double)gem_profile.profile[gtid].counters[counter]) : 0.0;
}
/*
 * Utils
 */
GEM_INLINE void PROF_SUM_REDUCE() {
  uint64_t i;
  for (i=1;i<gem_profile.num_threads;++i) {
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      gem_profile.profile[0].timers[j] += gem_profile.profile[i].timers[j];
      gem_profile.profile[0].counters[j] += gem_profile.profile[i].counters[j];
    }
  }
}
GEM_INLINE void PROF_SUM_OVERLAP() {
  uint64_t i;
  for (i=1;i<gem_profile.num_threads;++i) {
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      gem_profile.profile[0].timers[j] =
          MAX(gem_profile.profile[0].timers[j],gem_profile.profile[i].timers[j]);
      gem_profile.profile[0].counters[j] += gem_profile.profile[i].counters[j];
    }
  }
}
#endif

