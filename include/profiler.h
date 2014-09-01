/*
 * PROJECT: GEMMapper
 * FILE: profiler.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple time/functional profiler module
 */

#ifndef PROFILER_H_
#define PROFILER_H_

#include "commons.h"
#include "mm.h"

#include "gcounters.h"

/*
 * Precise Counter (enables calculation of mean,var,stddev at the cost of a little overhead)
 */
#define GEM_PRECISE_COUNTER

/*
 * Profiler Printers
 */
#define PRIcounter "lu(#%lu,m%lu,M%lu,{%lu})"
#define PRIcounter_content(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter),COUNTER_GET_MEAN(counter)
#define PRIcounterFull "lu(#%lu,m%lu,M%lu,{%lu,%lu,%lu})"
#define PRIcounterFull_content(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter), \
  COUNTER_GET_MEAN(counter),COUNTER_GET_VARIANCE(counter),COUNTER_GET_STDDEV(counter)

/*
 * Counters (from http://www.johndcook.com/standard_deviation.html)
 */
typedef struct {
  uint64_t total;
  uint64_t samples;
#ifdef GEM_PRECISE_COUNTER
  uint64_t min;
  uint64_t max;
  double m_oldM;
  double m_newM;
  double m_oldS;
  double m_newS;
#endif
} gem_counter_t;

GEM_INLINE void COUNTER_RESET(gem_counter_t* const counter);
GEM_INLINE void COUNTER_ADD(gem_counter_t* const counter,const uint64_t amount);

GEM_INLINE uint64_t COUNTER_GET_TOTAL(gem_counter_t* const counter);
GEM_INLINE uint64_t COUNTER_GET_NUM_SAMPLES(gem_counter_t* const counter);

GEM_INLINE uint64_t COUNTER_GET_MIN(gem_counter_t* const counter);
GEM_INLINE uint64_t COUNTER_GET_MAX(gem_counter_t* const counter);

GEM_INLINE double COUNTER_GET_MEAN(gem_counter_t* const counter);
GEM_INLINE double COUNTER_GET_VARIANCE(gem_counter_t* const counter);
GEM_INLINE double COUNTER_GET_STDDEV(gem_counter_t* const counter);

GEM_INLINE void COUNTER_COMBINE(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);

/*
 * Timers
 */
typedef struct {
  /* Timer */
  struct timespec begin_timer; // Timer begin
  struct timespec end_timer;   // Timer end
  /* Total time & samples taken */
  gem_counter_t time_ns;
  uint64_t accumulated;
} gem_timer_t;

GEM_INLINE void TIMER_START(gem_timer_t* const timer);
GEM_INLINE void TIMER_STOP(gem_timer_t* const timer);
GEM_INLINE void TIMER_PAUSE(gem_timer_t* const timer);
GEM_INLINE void TIMER_CONTINUE(gem_timer_t* const timer);
GEM_INLINE void TIMER_RESET(gem_timer_t* const timer);
GEM_INLINE void TIMER_RESTART(gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_NS(gem_timer_t* const timer);
GEM_INLINE double TIMER_GET_S(gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_MEAN(gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_VARIANCE(gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_STDDEV(gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_NUM_SAMPLES(gem_timer_t* const timer);

/*
 * Profiler
 */
#ifndef GEM_NOPROFILE /* GEM_PROFILE ENABLED */
// Constants
#define GP_MAX_COUNTERS 1000
// Profile Block
#define PROF_BLOCK()


/*
 * Setup
 */
GEM_INLINE void PROF_NEW(const uint64_t num_threads);
GEM_INLINE void PROF_DELETE();

/*
 * PROFILE-TIME functions
 */
GEM_INLINE void PROF_START_TIMER(const uint64_t timer);
GEM_INLINE void PROF_STOP_TIMER(const uint64_t timer);
GEM_INLINE void PROF_PAUSE_TIMER(const uint64_t timer);
GEM_INLINE void PROF_CONTINUE_TIMER(const uint64_t timer);
GEM_INLINE void PROF_RESET_TIMER(const uint64_t timer);
GEM_INLINE double PROF_GET_TIMER(const uint64_t timer);

/*
 * PROFILE-COUNTERS functions
 */
GEM_INLINE void PROF_RESET_COUNTER(const uint64_t counter);
GEM_INLINE uint64_t PROF_GET_COUNTER(const uint64_t counter);
GEM_INLINE void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value);
GEM_INLINE void PROF_INC_COUNTER(const uint64_t counter);

/*
 * Display statistics
 */
GEM_INLINE uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter);
GEM_INLINE double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2);
GEM_INLINE double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t tota_timer);
GEM_INLINE double PROF_TIME_PER_CALL(const uint64_t timer);

/*
 * Utils
 */
GEM_INLINE void PROF_SUM_REDUCE();

#else /* GEM_PROFILE DISABLED */
// Profile Block
#define PROF_BLOCK() if (0)
// Setup
#define PROF_NEW(num_threads);
#define PROF_DELETE();
// TIME functions
#define PROF_START_TIMER(timer);
#define PROF_STOP_TIMER(timer);
#define PROF_RESET_TIMER(timer);
#define PROF_GET_TIMER(timer);
// GENERAL COUNTERS functions
#define PROF_RESET_COUNTER(counter);
#define PROF_GET_COUNTER(counter);
#define PROF_SET_COUNTER(counter,value);
#define PROF_ADD_COUNTER(counter,value);
#define PROF_INC_COUNTER(counter);
#define PROF_DEC_COUNTER(counter);
// AGGREGATED functions
#define PROF_BEGIN(prof_counter);
#define PROF_END(prof_counter);
// Display statistics
#define PROF_TIME_PERCENTAGE(timer,tota_timer);
#define PROF_COUNT_PERCENTAGE(counter,total_counter);
#define PROF_COUNT_DIV(counter1,counter2);
#define PROF_TIME_PER_CALL(timer,counter);
// Utils
#define PROF_SUM_REDUCE();
#define PROF_SUM_OVERLAP();
#endif /* GEM_NOPROFILE */

#endif /* PROFILER_H_ */
