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

/*
 * Precise Counter (enables calculation of mean,var,stddev at the cost of a little overhead)
 */
#define GEM_PRECISE_COUNTER

/*
 * Profiler Printers
 */
// Counters
#define PRIcounter "lu(#%lu,m%lu,M%lu,{%.2f})"
#define PRIcounterVal(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter),COUNTER_GET_MEAN(counter)
#ifdef GEM_PRECISE_COUNTER
#define PRIcounterX "lu(#%lu,m%lu,M%lu,{%.2f,%.2f,%.2f})"
#define PRIcounterXVal(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter), \
  COUNTER_GET_MEAN(counter),COUNTER_GET_VARIANCE(counter),COUNTER_GET_STDDEV(counter)
#else
#define PRIcounterX             PRIcounter
#define PRIcounterXVal(counter) PRIcounterVal(counter)
#endif

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

GEM_INLINE uint64_t COUNTER_GET_TOTAL(const gem_counter_t* const counter);
GEM_INLINE uint64_t COUNTER_GET_NUM_SAMPLES(const gem_counter_t* const counter);
GEM_INLINE uint64_t COUNTER_GET_MIN(const gem_counter_t* const counter);
GEM_INLINE uint64_t COUNTER_GET_MAX(const gem_counter_t* const counter);
GEM_INLINE double COUNTER_GET_MEAN(const gem_counter_t* const counter);
GEM_INLINE double COUNTER_GET_VARIANCE(const gem_counter_t* const counter);
GEM_INLINE double COUNTER_GET_STDDEV(const gem_counter_t* const counter);

GEM_INLINE void COUNTER_COMBINE(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);

GEM_INLINE void COUNTER_PRINT(
    FILE* const stream,const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,const char* const units,const bool full_report);
GEM_INLINE void PERCENTAGE_PRINT(FILE* const stream,const gem_counter_t* const counter);

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

#define TIMER_CONVERT_NS_TO_US(time_ns) ((double)(time_ns)/1E3)
#define TIMER_CONVERT_NS_TO_MS(time_ns) ((double)(time_ns)/1E6)
#define TIMER_CONVERT_NS_TO_S(time_ns)  ((double)(time_ns)/1E9)
#define TIMER_CONVERT_NS_TO_M(time_ns)  ((double)(time_ns)/1E9/60.0)
#define TIMER_CONVERT_NS_TO_H(time_ns)  ((double)(time_ns)/1E9/3600.0)

GEM_INLINE void TIMER_START(gem_timer_t* const timer);
GEM_INLINE void TIMER_STOP(gem_timer_t* const timer);
GEM_INLINE void TIMER_PAUSE(gem_timer_t* const timer);
GEM_INLINE void TIMER_CONTINUE(gem_timer_t* const timer);
GEM_INLINE void TIMER_RESET(gem_timer_t* const timer);
GEM_INLINE void TIMER_RESTART(gem_timer_t* const timer);

GEM_INLINE uint64_t TIMER_GET_TOTAL_NS(const gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_NUM_SAMPLES(const gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_MIN_NS(const gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_MAX_NS(const gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_MEAN(const gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_VARIANCE(const gem_timer_t* const timer);
GEM_INLINE uint64_t TIMER_GET_STDDEV(const gem_timer_t* const timer);

GEM_INLINE void TIMER_PRINT(
    FILE* const stream,const gem_timer_t* const timer,const gem_timer_t* const ref_timer);

#define TIMER_GET_TOTAL_US(timer) TIMER_CONVERT_NS_TO_US(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_MS(timer) TIMER_CONVERT_NS_TO_MS(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_S(timer)  TIMER_CONVERT_NS_TO_S(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_M(timer)  TIMER_CONVERT_NS_TO_M(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_H(timer)  TIMER_CONVERT_NS_TO_H(TIMER_GET_TOTAL_NS(timer))

/*
 * Reference Counter (Counts wrt a reference counter. Eg ranks)
 */
typedef struct {
  /* Rank counter */
  uint64_t begin_count;
  /* Total ranks & samples taken */
  gem_counter_t counter;
  uint64_t accumulated;
} gem_reference_counter_t;

GEM_INLINE void RCOUNTER_START(gem_reference_counter_t* const rcounter,const uint64_t reference);
GEM_INLINE void RCOUNTER_STOP(gem_reference_counter_t* const rcounter,const uint64_t reference);
GEM_INLINE void RCOUNTER_PAUSE(gem_reference_counter_t* const rcounter,const uint64_t reference);
GEM_INLINE void RCOUNTER_CONTINUE(gem_reference_counter_t* const rcounter,const uint64_t reference);
GEM_INLINE void RCOUNTER_RESET(gem_reference_counter_t* const rcounter);
GEM_INLINE void RCOUNTER_RESTART(gem_reference_counter_t* const rcounter);
GEM_INLINE uint64_t RCOUNTER_GET_TOTAL(gem_reference_counter_t* const rcounter);
GEM_INLINE uint64_t RCOUNTER_GET_NUM_SAMPLES(gem_reference_counter_t* const rcounter);
GEM_INLINE uint64_t RCOUNTER_GET_MIN(gem_reference_counter_t* const rcounter);
GEM_INLINE uint64_t RCOUNTER_GET_MAX(gem_reference_counter_t* const rcounter);
GEM_INLINE uint64_t RCOUNTER_GET_MEAN(gem_reference_counter_t* const rcounter);
GEM_INLINE uint64_t RCOUNTER_GET_VARIANCE(gem_reference_counter_t* const rcounter);
GEM_INLINE uint64_t RCOUNTER_GET_STDDEV(gem_reference_counter_t* const rcounter);

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
GEM_INLINE gem_timer_t* PROF_GET_TIMER(const uint64_t timer);

/*
 * PROFILE-COUNTERS functions
 */
GEM_INLINE void PROF_RESET_COUNTER(const uint64_t counter);
GEM_INLINE void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value);
GEM_INLINE void PROF_INC_COUNTER(const uint64_t counter);
GEM_INLINE gem_counter_t* PROF_GET_COUNTER(const uint64_t counter);

/*
 * PROFILE-RANKS functions
 */
GEM_INLINE void PROF_START_RANK(const uint64_t rank);
GEM_INLINE void PROF_STOP_RANK(const uint64_t rank);
GEM_INLINE void PROF_PAUSE_RANK(const uint64_t rank);
GEM_INLINE void PROF_CONTINUE_RANK(const uint64_t rank);
GEM_INLINE void PROF_RESET_RANK(const uint64_t rank);
GEM_INLINE gem_counter_t* PROF_GET_RANK(const uint64_t rank);

/*
 * PROFILE-COMBINED (TIME/RANKS) functions
 */
GEM_INLINE void PROF_START(const uint64_t timer__ranks);
GEM_INLINE void PROF_STOP(const uint64_t timer__ranks);
GEM_INLINE void PROF_PAUSE(const uint64_t timer__ranks);
GEM_INLINE void PROF_CONTINUE(const uint64_t timer__ranks);
GEM_INLINE void PROF_RESET(const uint64_t timer__ranks);

/*
 * Display statistics
 */
GEM_INLINE uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter);
GEM_INLINE double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2);
GEM_INLINE double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t total_timer);
GEM_INLINE double PROF_TIME_PER_CALL(const uint64_t timer); // ms/call

/*
 * Utils
 */
GEM_INLINE void PROF_SUM_REDUCE();

#else /* GEM_PROFILE DISABLED */
// Profile Block
#define PROF_BLOCK() if (0)
// TODO
#endif /* GEM_NOPROFILE */

#endif /* PROFILER_H_ */
