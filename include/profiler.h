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
#include "profiler_cuda.h"

/*
 * Profiler Printers
 */
// Counters
#define PRIcounter "lu(#%"PRIu64",m%"PRIu64",M%"PRIu64",{%.2f})"
#define PRIcounterVal(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter),COUNTER_GET_MEAN(counter)
#define PRIcounterX "lu(#%"PRIu64",m%"PRIu64",M%"PRIu64",{%.2f,%.2f,%.2f})"
#define PRIcounterXVal(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter), \
  COUNTER_GET_MEAN(counter),COUNTER_GET_VARIANCE(counter),COUNTER_GET_STDDEV(counter)

/*
 * Counters (from http://www.johndcook.com/standard_deviation.html)
 */
typedef struct {
  uint64_t total;
  uint64_t samples;
  uint64_t min;
  uint64_t max;
  double m_oldM;
  double m_newM;
  double m_oldS;
  double m_newS;
} gem_counter_t;

void COUNTER_RESET(gem_counter_t* const counter);
void COUNTER_ADD(gem_counter_t* const counter,const uint64_t amount);

uint64_t COUNTER_GET_TOTAL(const gem_counter_t* const counter);
uint64_t COUNTER_GET_NUM_SAMPLES(const gem_counter_t* const counter);
uint64_t COUNTER_GET_MIN(const gem_counter_t* const counter);
uint64_t COUNTER_GET_MAX(const gem_counter_t* const counter);
double COUNTER_GET_MEAN(const gem_counter_t* const counter);
double COUNTER_GET_VARIANCE(const gem_counter_t* const counter);
double COUNTER_GET_STDDEV(const gem_counter_t* const counter);

void COUNTER_COMBINE_MAX(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);
void COUNTER_COMBINE_MIN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);
void COUNTER_COMBINE_MEAN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);

void COUNTER_PRINT(
    FILE* const stream,const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,const char* const units,const bool full_report);
void SAMPLER_PRINT(
    FILE* const stream,const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,const char* const units);
void PERCENTAGE_PRINT(FILE* const stream,const gem_counter_t* const counter);

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

void TIMER_START(gem_timer_t* const timer);
void TIMER_STOP(gem_timer_t* const timer);
void TIMER_PAUSE(gem_timer_t* const timer);
void TIMER_CONTINUE(gem_timer_t* const timer);
void TIMER_RESET(gem_timer_t* const timer);
void TIMER_RESTART(gem_timer_t* const timer);

uint64_t TIMER_GET_TOTAL_NS(const gem_timer_t* const timer);
uint64_t TIMER_GET_NUM_SAMPLES(const gem_timer_t* const timer);
uint64_t TIMER_GET_MIN_NS(const gem_timer_t* const timer);
uint64_t TIMER_GET_MAX_NS(const gem_timer_t* const timer);
uint64_t TIMER_GET_MEAN(const gem_timer_t* const timer);
uint64_t TIMER_GET_VARIANCE(const gem_timer_t* const timer);
uint64_t TIMER_GET_STDDEV(const gem_timer_t* const timer);

void TIMER_PRINT(
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

void RCOUNTER_START(gem_reference_counter_t* const rcounter,const uint64_t reference);
void RCOUNTER_STOP(gem_reference_counter_t* const rcounter,const uint64_t reference);
void RCOUNTER_PAUSE(gem_reference_counter_t* const rcounter,const uint64_t reference);
void RCOUNTER_CONTINUE(gem_reference_counter_t* const rcounter,const uint64_t reference);
void RCOUNTER_RESET(gem_reference_counter_t* const rcounter);
void RCOUNTER_RESTART(gem_reference_counter_t* const rcounter);
uint64_t RCOUNTER_GET_TOTAL(gem_reference_counter_t* const rcounter);
uint64_t RCOUNTER_GET_NUM_SAMPLES(gem_reference_counter_t* const rcounter);
uint64_t RCOUNTER_GET_MIN(gem_reference_counter_t* const rcounter);
uint64_t RCOUNTER_GET_MAX(gem_reference_counter_t* const rcounter);
uint64_t RCOUNTER_GET_MEAN(gem_reference_counter_t* const rcounter);
uint64_t RCOUNTER_GET_VARIANCE(gem_reference_counter_t* const rcounter);
uint64_t RCOUNTER_GET_STDDEV(gem_reference_counter_t* const rcounter);

/*
 * Profiler
 */
// #define GEM_NOPROFILE /* Disables profile */
#define GP_MAX_COUNTERS 1000
typedef enum { reduce_sum, reduce_max, reduce_min, reduce_mean, reduce_sample } profile_reduce_type;

#ifndef GEM_NOPROFILE /* GEM_PROFILE ENABLED */

/*
 * Profile Block
 */
#define PROF_BLOCK()

/*
 * Setup
 */
void PROF_NEW(const uint64_t num_threads);
void PROF_DELETE();

/*
 * PROFILE-TIME functions
 */
void PROF_START_TIMER(const uint64_t timer);
void PROF_STOP_TIMER(const uint64_t timer);
void PROF_PAUSE_TIMER(const uint64_t timer);
void PROF_CONTINUE_TIMER(const uint64_t timer);
void PROF_RESET_TIMER(const uint64_t timer);
gem_timer_t* PROF_GET_TIMER(const uint64_t timer);

/*
 * PROFILE-COUNTERS functions
 */
void PROF_RESET_COUNTER(const uint64_t counter);
void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value);
void PROF_INC_COUNTER(const uint64_t counter);
gem_counter_t* PROF_GET_COUNTER(const uint64_t counter);

/*
 * PROFILE-RANKS functions
 */
void PROF_START_RANK(const uint64_t rank);
void PROF_STOP_RANK(const uint64_t rank);
void PROF_PAUSE_RANK(const uint64_t rank);
void PROF_CONTINUE_RANK(const uint64_t rank);
void PROF_RESET_RANK(const uint64_t rank);
gem_counter_t* PROF_GET_RANK(const uint64_t rank);

/*
 * PROFILE-COMBINED (TIME/RANKS) functions
 */
#define PROF_START(TAG) \
  PROFILE_CUDA_START(#TAG,TAG%6); \
  _PROF_START(TAG)
#define PROF_STOP(TAG) \
  PROFILE_CUDA_STOP; \
  _PROF_STOP(TAG)
void _PROF_START(const uint64_t timer__ranks);
void _PROF_STOP(const uint64_t timer__ranks);
void PROF_PAUSE(const uint64_t timer__ranks);
void PROF_CONTINUE(const uint64_t timer__ranks);
void PROF_RESET(const uint64_t timer__ranks);

/*
 * Display statistics
 */
uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter);
double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2);
double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t total_timer);
double PROF_TIME_PER_CALL(const uint64_t timer); // ms/call

/*
 * Utils
 */
void PROF_REDUCE_SUM();
void PROF_REDUCE_MAX();
void PROF_REDUCE_MIN();
void PROF_REDUCE_MEAN();
void PROF_REDUCE_SAMPLE();

#else /* GEM_PROFILE DISABLED */
// Profile Block
#define PROF_BLOCK() if (0)
// Setup
#define PROF_NEW(num_threads);
#define PROF_DELETE();
// PROFILE-TIME functions
#define PROF_START_TIMER(timer);
#define PROF_STOP_TIMER(timer);
#define PROF_PAUSE_TIMER(timer);
#define PROF_CONTINUE_TIMER(timer);
#define PROF_RESET_TIMER(timer);
#define PROF_GET_TIMER(timer);
// PROFILE-COUNTERS functions
#define PROF_RESET_COUNTER(counter);
#define PROF_ADD_COUNTER(counter,value);
#define PROF_INC_COUNTER(counter);
#define PROF_GET_COUNTER(counter);
// PROFILE-RANKS functions
#define PROF_START_RANK(rank);
#define PROF_STOP_RANK(rank);
#define PROF_PAUSE_RANK(rank);
#define PROF_CONTINUE_RANK(rank);
#define PROF_RESET_RANK(rank);
#define PROF_GET_RANK(rank);
// PROFILE-COMBINED (TIME/RANKS) functions
#define PROF_START(timer__ranks);
#define PROF_STOP(timer__ranks);
#define PROF_PAUSE(timer__ranks);
#define PROF_CONTINUE(timer__ranks);
#define PROF_RESET(timer__ranks);
// Display statistics
#define PROF_COUNT_PERCENTAGE(counter,total_counter)
#define PROF_COUNT_DIV(counter1,counter2)
#define PROF_TIME_PERCENTAGE(timer,total_timer)
#define PROF_TIME_PER_CALL(timer)
// Utils
#define PROF_REDUCE_SUM();
#define PROF_REDUCE_MAX();
#define PROF_REDUCE_MIN();
#define PROF_REDUCE_MEAN();
#define PROF_REDUCE_SAMPLE();
#endif /* GEM_NOPROFILE */

#endif /* PROFILER_H_ */
