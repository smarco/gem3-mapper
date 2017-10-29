/*
 * PROJECT: GEM-Tools library
 * FILE: gt_profile.c
 * DATE: 16/03/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_PROFILER_H_
#define GT_PROFILER_H_

#include "gt_commons.h"
#include "gt_mm.h"

/*
 * Profiler
 */
typedef struct {
  // Time counters
  struct timeval* gt_prof_begin_timer;
  struct timeval* gt_prof_end_timer;
  double* gt_prof_timer;
  // General counters & functions
  uint64_t *gt_prof_counter;
  // Maximum number of counters available
  uint64_t max_num_counters;
} gt_profile;

/*
 * Setup
 */
GT_INLINE gt_profile* GPROF_NEW(const uint64_t max_size);
GT_INLINE void GPROF_DELETE(gt_profile* const profile);

/*
 * TIME functions
 */
GT_INLINE void GPROF_START_TIMER(gt_profile* const profile,const uint64_t timer);
GT_INLINE void GPROF_STOP_TIMER(gt_profile* const profile,const uint64_t timer);
GT_INLINE void GPROF_RESET_TIMER(gt_profile* const profile,const uint64_t timer);
GT_INLINE double GPROF_GET_TIMER(gt_profile* const profile,const uint64_t timer);

/*
 * GENERAL COUNTERS functions
 */
GT_INLINE void GPROF_RESET_COUNTER(gt_profile* const profile,const uint64_t counter);

GT_INLINE void GPROF_SET_COUNTER(gt_profile* const profile,const uint64_t counter,const uint64_t value);
GT_INLINE void GPROF_ADD_COUNTER(gt_profile* const profile,const uint64_t counter,const uint64_t value);
GT_INLINE void GPROF_INC_COUNTER(gt_profile* const profile,const uint64_t counter);
GT_INLINE void GPROF_DEC_COUNTER(gt_profile* const profile,const uint64_t counter);

GT_INLINE uint64_t GPROF_GET_COUNTER(gt_profile* const profile,const uint64_t counter);

/*
 * Display statistics
 */
GT_INLINE double GPROF_TIME_PERCENTAGE(gt_profile* const profile,const uint64_t timer,const uint64_t tota_timer);
GT_INLINE uint64_t GPROF_COUNT_PERCENTAGE(gt_profile* const profile,const uint64_t counter,const uint64_t total_counter);
GT_INLINE double GPROF_COUNT_DIV(gt_profile* const profile,const uint64_t counter1,const uint64_t counter2);
GT_INLINE double GPROF_TIME_PER_CALL(gt_profile* const profile,const uint64_t timer,const uint64_t counter);

/*
 * Utils
 */
GT_INLINE void GPROF_SUM_REDUCE(gt_profile** const profile,const uint64_t num_profiles);
GT_INLINE void GPROF_SUM_OVERLAP(gt_profile** const profile,const uint64_t num_profiles);

#endif /* GT_PROFILER_H_ */
