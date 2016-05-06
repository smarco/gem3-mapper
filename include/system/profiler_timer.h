/*
 * PROJECT: GEMMapper
 * FILE: profiler_timer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#ifndef PROFILER_TIMER_H
#define PROFILER_TIMER_H

#include "system/commons.h"
#include "stats/stats_counter.h"

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
    FILE* const stream,const gem_timer_t* const timer,
    const gem_timer_t* const ref_timer);

#define TIMER_GET_TOTAL_US(timer) TIMER_CONVERT_NS_TO_US(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_MS(timer) TIMER_CONVERT_NS_TO_MS(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_S(timer)  TIMER_CONVERT_NS_TO_S(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_M(timer)  TIMER_CONVERT_NS_TO_M(TIMER_GET_TOTAL_NS(timer))
#define TIMER_GET_TOTAL_H(timer)  TIMER_CONVERT_NS_TO_H(TIMER_GET_TOTAL_NS(timer))

#endif /* PROFILER_TIMER_H */
