/*
 * PROJECT: GEM-Tools library
 * FILE: gt_profile.c
 * DATE: 16/03/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_profiler.h"

// Operative macros
#define GT_TIME_DIFF(start,end) ((end.tv_sec + end.tv_usec/1E6) - (start.tv_sec + start.tv_usec/1E6))

/*
 * Setup
 */
GT_INLINE gt_profile* GPROF_NEW(const uint64_t max_num_counters) {
  gt_profile* const profile = gt_alloc(gt_profile);
  // Allocate all profile counters
  profile->gt_prof_begin_timer = gt_calloc(max_num_counters,struct timeval,true);
  profile->gt_prof_end_timer = gt_calloc(max_num_counters,struct timeval,true);
  profile->gt_prof_timer = gt_calloc(max_num_counters,double,true);
  profile->gt_prof_counter = gt_calloc(max_num_counters,uint64_t,true);
  profile->max_num_counters = max_num_counters;
  return profile;
}
GT_INLINE void GPROF_DELETE(gt_profile* const profile) {
  // Release all profile counters
  gt_free(profile->gt_prof_begin_timer);
  gt_free(profile->gt_prof_end_timer);
  gt_free(profile->gt_prof_timer);
  gt_free(profile->gt_prof_counter);
}
/*
 * TIME functions
 */
GT_INLINE void GPROF_START_TIMER(gt_profile* const profile,const uint64_t timer) {
  gettimeofday((profile->gt_prof_begin_timer+timer),NULL);
}
GT_INLINE void GPROF_STOP_TIMER(gt_profile* const profile,const uint64_t timer) {
  gettimeofday((profile->gt_prof_end_timer+timer),NULL);
  profile->gt_prof_timer[timer]+=GT_TIME_DIFF(profile->gt_prof_begin_timer[timer],profile->gt_prof_end_timer[timer]);
}
GT_INLINE void GPROF_RESET_TIMER(gt_profile* const profile,const uint64_t timer) {
  profile->gt_prof_timer[timer]=0.0;
}
GT_INLINE double GPROF_GET_TIMER(gt_profile* const profile,const uint64_t timer) {
  return profile->gt_prof_timer[timer];
}
/*
 * GENERAL COUNTERS functions
 */
GT_INLINE void GPROF_RESET_COUNTER(gt_profile* const profile,const uint64_t counter) {
  profile->gt_prof_counter[counter] = 0;
}
GT_INLINE void GPROF_SET_COUNTER(gt_profile* const profile,const uint64_t counter,const uint64_t value) {
  profile->gt_prof_counter[counter] = value;
}
GT_INLINE void GPROF_ADD_COUNTER(gt_profile* const profile,const uint64_t counter,const uint64_t value) {
  profile->gt_prof_counter[counter] += value;
}
GT_INLINE void GPROF_INC_COUNTER(gt_profile* const profile,const uint64_t counter) {
  ++(profile->gt_prof_counter[counter]);
}
GT_INLINE void GPROF_DEC_COUNTER(gt_profile* const profile,const uint64_t counter) {
  --(profile->gt_prof_counter[counter]);
}
GT_INLINE uint64_t GPROF_GET_COUNTER(gt_profile* const profile,const uint64_t counter) {
  return profile->gt_prof_counter[counter];
}
/*
 * Display statistics
 */
GT_INLINE double GPROF_TIME_PERCENTAGE(gt_profile* const profile,const uint64_t timer,const uint64_t tota_timer) {
  return GT_GET_PERCENTAGE(profile->gt_prof_timer[timer],profile->gt_prof_timer[tota_timer]);
}
GT_INLINE uint64_t GPROF_COUNT_PERCENTAGE(gt_profile* const profile,const uint64_t counter,const uint64_t total_counter) {
  return GT_GET_PERCENTAGE(profile->gt_prof_counter[counter],profile->gt_prof_counter[total_counter]);
}
GT_INLINE double GPROF_COUNT_DIV(gt_profile* const profile,const uint64_t counter1,const uint64_t counter2) {
  return (float)profile->gt_prof_counter[counter1]/(float)profile->gt_prof_counter[counter2];
}
GT_INLINE double GPROF_TIME_PER_CALL(gt_profile* const profile,const uint64_t timer,const uint64_t counter) {
  return (profile->gt_prof_counter[counter]!=0) ?
      1000.0*((double)profile->gt_prof_timer[timer]/(double)profile->gt_prof_counter[counter]) : 0.0;
}
/*
 * Utils
 */
GT_INLINE void GPROF_SUM_REDUCE(gt_profile** const profile,const uint64_t num_profiles) {
  uint64_t i;
  for (i=1;i<num_profiles;++i) {
    uint64_t j;
    for (j=0;j<profile[0]->max_num_counters;++j) {
      profile[0]->gt_prof_timer[j] += profile[i]->gt_prof_timer[j];
      profile[0]->gt_prof_counter[j] += profile[i]->gt_prof_counter[j];
    }
  }
}
GT_INLINE void GPROF_SUM_OVERLAP(gt_profile** const profile,const uint64_t num_profiles) {
  uint64_t i;
  for (i=1;i<num_profiles;++i) {
    uint64_t j;
    for (j=0;j<profile[0]->max_num_counters;++j) {
      profile[0]->gt_prof_timer[j] = GT_MAX(profile[0]->gt_prof_timer[j],profile[i]->gt_prof_timer[j]);
      profile[0]->gt_prof_counter[j] += profile[i]->gt_prof_counter[j];
    }
  }
}


