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
 * Timers
 */
typedef struct {
  uint64_t time_acc;    // Time accumulated
  struct timeval begin_timer; // Timer begin
  struct timeval end_timer;   // Timer end
} gem_timer_t;

GEM_INLINE void START_TIMER(gem_timer_t* const timer);
GEM_INLINE void STOP_TIMER(gem_timer_t* const timer);
GEM_INLINE void RESET_TIMER(gem_timer_t* const timer);
GEM_INLINE void RESET__START_TIMER(gem_timer_t* const timer);
GEM_INLINE double GET_TIMER(gem_timer_t* const timer);

#ifndef GEM_NOPROFILE /* GEM_PROFILE ENABLED */
/*
 * Setup
 */
GEM_INLINE void PROF_NEW(const uint64_t num_threads);
GEM_INLINE void PROF_DELETE();

/*
 * TIME functions
 */
GEM_INLINE void PROF_START_TIMER(const uint64_t timer);
GEM_INLINE void PROF_STOP_TIMER(const uint64_t timer);
GEM_INLINE void PROF_RESET_TIMER(const uint64_t timer);
GEM_INLINE double PROF_GET_TIMER(const uint64_t timer);

/*
 * GENERAL COUNTERS functions
 */
GEM_INLINE void PROF_RESET_COUNTER(const uint64_t counter);

GEM_INLINE uint64_t PROF_GET_COUNTER(const uint64_t counter);
GEM_INLINE void PROF_SET_COUNTER(const uint64_t counter,const uint64_t value);

GEM_INLINE void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value);
GEM_INLINE void PROF_INC_COUNTER(const uint64_t counter);
GEM_INLINE void PROF_DEC_COUNTER(const uint64_t counter);

/*
 * AGGREGATED functions
 */
GEM_INLINE void PROF_BEGIN(const uint64_t prof_counter);
GEM_INLINE void PROF_END(const uint64_t prof_counter);

/*
 * Display statistics
 */
GEM_INLINE double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t tota_timer);
GEM_INLINE uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter);
GEM_INLINE double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2);
GEM_INLINE double PROF_TIME_PER_CALL(const uint64_t timer,const uint64_t counter);

/*
 * Utils
 */
GEM_INLINE void PROF_SUM_REDUCE();
GEM_INLINE void PROF_SUM_OVERLAP();

#else /* GEM_PROFILE DISABLED */
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
