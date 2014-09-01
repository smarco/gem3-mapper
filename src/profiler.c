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
 * Counters
 */

GEM_INLINE void COUNTER_RESET(gem_counter_t* const counter) {
  counter->total = 0.0;
  counter->samples = 0;
}
GEM_INLINE void COUNTER_ADD(gem_counter_t* const counter,const uint64_t amount) {
  // Add to total & increment number of samples
  counter->total += amount;
  ++(counter->samples);
#ifdef GEM_PRECISE_COUNTER
  counter->min = MIN(counter->min,amount);
  counter->max = MAX(counter->max,amount);
  // See Knuth TAOCP vol 2, 3rd edition, page 232
  if (counter->samples == 1) {
    counter->m_oldM = amount;
    counter->m_newM = amount;
    counter->m_oldS = 0.0;
  } else {
    counter->m_newM = counter->m_oldM + ((double)amount-counter->m_oldM)/(double)counter->samples;
    counter->m_newS = counter->m_oldS + ((double)amount-counter->m_oldM)*((double)amount-counter->m_newM);
    counter->m_oldM = counter->m_newM;
    counter->m_oldS = counter->m_newS;
  }
#endif
}
GEM_INLINE uint64_t COUNTER_GET_TOTAL(gem_counter_t* const counter) {
  return counter->total;
}
GEM_INLINE uint64_t COUNTER_GET_NUM_SAMPLES(gem_counter_t* const counter) {
  return counter->samples;
}
GEM_INLINE uint64_t COUNTER_GET_MIN(gem_counter_t* const counter) {
#ifdef GEM_PRECISE_COUNTER
  return counter->min;
#else
  return 0;
#endif
}
GEM_INLINE uint64_t COUNTER_GET_MAX(gem_counter_t* const counter) {
#ifdef GEM_PRECISE_COUNTER
  return counter->max;
#else
  return 0;
#endif
}
GEM_INLINE double COUNTER_GET_MEAN(gem_counter_t* const counter) {
#ifdef GEM_PRECISE_COUNTER
  return (counter->samples > 0) ? counter->m_newM : 0.0;
#else
  return (double)counter->total/(double)counter->samples;
#endif
}
GEM_INLINE double COUNTER_GET_VARIANCE(gem_counter_t* const counter) {
#ifdef GEM_PRECISE_COUNTER
  return ((counter->samples > 1) ? counter->m_newS/(double)(counter->samples - 1) : 0.0);
#else
  return 0.0;
#endif
}
GEM_INLINE double COUNTER_GET_STDDEV(gem_counter_t* const counter) {
#ifdef GEM_PRECISE_COUNTER
  return sqrt(COUNTER_GET_VARIANCE(counter));
#else
  return 0.0;
#endif
}
GEM_INLINE void COUNTER_COMBINE(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
  counter_dst->total += counter_src->total;
  counter_dst->samples += counter_src->samples;
#ifdef GEM_PRECISE_COUNTER
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
  counter_dst->m_newS = 0.0; // TODO
  counter_dst->m_newM = 0.0;
  counter_dst->m_oldS = 0.0;
  counter_dst->m_oldM = 0.0;
#endif
}
/*
 * Timers
 */
GEM_INLINE void TIMER_START(gem_timer_t* const timer) {
  timer->accumulated = 0;
  TIMER_CONTINUE(timer);
}
GEM_INLINE void TIMER_STOP(gem_timer_t* const timer) {
  TIMER_PAUSE(timer);
  COUNTER_ADD(&timer->time_ns,timer->accumulated);
}
GEM_INLINE void TIMER_PAUSE(gem_timer_t* const timer) {
  clock_gettime(CLOCK_REALTIME,&timer->end_timer);
  timer->accumulated += TIME_DIFF_NS(timer->begin_timer,timer->end_timer);
}
GEM_INLINE void TIMER_CONTINUE(gem_timer_t* const timer) {
  clock_gettime(CLOCK_REALTIME,&timer->begin_timer);
}
GEM_INLINE void TIMER_RESET(gem_timer_t* const timer) {
  COUNTER_RESET(&timer->time_ns);
}
GEM_INLINE void TIMER_RESTART(gem_timer_t* const timer) {
  TIMER_RESET(timer);
  TIMER_START(timer);
}
GEM_INLINE uint64_t TIMER_GET_NS(gem_timer_t* const timer) {
  return COUNTER_GET_TOTAL(&timer->time_ns);
}
GEM_INLINE double TIMER_GET_S(gem_timer_t* const timer) {
  return (double)COUNTER_GET_TOTAL(&timer->time_ns)/1E9;
}
GEM_INLINE uint64_t TIMER_GET_MEAN(gem_timer_t* const timer) {
  return COUNTER_GET_MEAN(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_VARIANCE(gem_timer_t* const timer) {
  return COUNTER_GET_VARIANCE(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_STDDEV(gem_timer_t* const timer) {
  return COUNTER_GET_STDDEV(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_NUM_SAMPLES(gem_timer_t* const timer) {
  return COUNTER_GET_NUM_SAMPLES(&timer->time_ns);
}

/*
 * Profile
 */
#ifndef GEM_NOPROFILE
typedef struct {
  gem_timer_t* timers; // Time counters
  gem_counter_t* counters; // General counters & functions
} profile_t;
typedef struct {
  // Profiler
  profile_t* profile;
  // Limits
  uint64_t num_threads;
} profiler_t;
profiler_t gem_profile; // THE GREAT PROFILER

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
    gem_profile.profile[i].timers = mm_calloc(GP_MAX_COUNTERS,gem_timer_t,true);
    gem_profile.profile[i].counters = mm_calloc(GP_MAX_COUNTERS,gem_counter_t,true);
  }
}
GEM_INLINE void PROF_DELETE() {
  // Release all profile counters
  uint64_t i;
  for (i=0;i<gem_profile.num_threads;++i) {
    mm_free(gem_profile.profile[i].timers);
    mm_free(gem_profile.profile[i].counters);
  }
  mm_free(gem_profile.profile);
}
/*
 * PROFILE-TIME functions
 */
GEM_INLINE void PROF_START_TIMER(const uint64_t timer) {
  TIMER_START(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_STOP_TIMER(const uint64_t timer) {
  TIMER_STOP(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_PAUSE_TIMER(const uint64_t timer) {
  TIMER_PAUSE(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_CONTINUE_TIMER(const uint64_t timer) {
  TIMER_CONTINUE(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_RESET_TIMER(const uint64_t timer) {
  TIMER_RESET(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE double PROF_GET_TIMER(const uint64_t timer) {
  return TIMER_GET_S(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE double PROF_GET_NUM_SAMPLES(const uint64_t timer) {
  return TIMER_GET_NUM_SAMPLES(gem_profile.profile[gtid()].timers+timer);
}
/*
 * PROFILE-COUNTERS functions
 */
GEM_INLINE void PROF_RESET_COUNTER(const uint64_t counter) {
  COUNTER_RESET(gem_profile.profile[gtid()].counters+counter);
}
GEM_INLINE uint64_t PROF_GET_COUNTER(const uint64_t counter) {
  return COUNTER_GET_TOTAL(gem_profile.profile[gtid()].counters+counter);
}
GEM_INLINE void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value) {
  COUNTER_ADD(gem_profile.profile[gtid()].counters+counter,value);
}
GEM_INLINE void PROF_INC_COUNTER(const uint64_t counter) {
  COUNTER_ADD(gem_profile.profile[gtid()].counters+counter,1);
}
/*
 * Display statistics
 */
GEM_INLINE uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter) {
  const uint64_t gtid = gtid();
  return PERCENTAGE(
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters + counter),
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters + total_counter) );
}
GEM_INLINE double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2) {
  const uint64_t gtid = gtid();
  return DIV(
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters+counter1),
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters+counter2) );
}
GEM_INLINE double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t tota_timer) {
  const uint64_t gtid = gtid();
  return PERCENTAGE(
      TIMER_GET_S(gem_profile.profile[gtid].timers + timer),
      TIMER_GET_S(gem_profile.profile[gtid].timers + tota_timer));
}
GEM_INLINE double PROF_TIME_PER_CALL(const uint64_t timer) {
  const uint64_t gtid = gtid();
  const uint64_t num_calls = TIMER_GET_NUM_SAMPLES(gem_profile.profile[gtid].timers+timer);
  return (num_calls!=0) ?
      1000.0*(TIMER_GET_S(gem_profile.profile[gtid].timers+timer) / num_calls) : 0.0;
}
/*
 * Utils
 */
GEM_INLINE void PROF_SUM_REDUCE() {
  uint64_t i;
  for (i=1;i<gem_profile.num_threads;++i) {
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      COUNTER_COMBINE(&gem_profile.profile[0].timers[j].time_ns,&gem_profile.profile[i].timers[j].time_ns);
      COUNTER_COMBINE(gem_profile.profile[0].counters+j,gem_profile.profile[i].counters+j);
    }
  }
}
#endif /* !GEM_NOPROFILE */

