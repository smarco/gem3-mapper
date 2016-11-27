/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "profiler/profiler_gem.h"
#include "system/gthread.h"

#ifdef GEM_PROFILE

/*
 * THE GREAT PROFILER
 */
profiler_t gem_profile;

/*
 * Setup
 */
void PROF_NEW(const uint64_t num_threads) {
  // Allocate handler
  gem_profile.profile = mm_calloc(num_threads,profile_t,true);
  // Initialize
  gem_profile.num_threads = num_threads;
  // Allocate all profiles
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    // Allocate
    gem_profile.profile[i].timers = mm_calloc(GP_MAX_COUNTERS,gem_timer_t,true);
    gem_profile.profile[i].counters = mm_calloc(GP_MAX_COUNTERS,gem_counter_t,true);
    gem_profile.profile[i].ranks = mm_calloc(GP_MAX_COUNTERS,gem_reference_counter_t,true);
    // Others
    gem_profile.profile[i].strata_deltas_edit = stats_vector_step_range_new(10,1,10);
    gem_profile.profile[i].strata_deltas_swg = stats_vector_step_range_new(20,1,10);
    // Initialize minimums (towards combine helps a lot)
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      gem_profile.profile[i].timers[j].time_ns.min = UINT64_MAX;
      gem_profile.profile[i].counters[j].min = UINT64_MAX;
      gem_profile.profile[i].ranks[j].counter.min = UINT64_MAX;
    }
  }
}
void PROF_DELETE(void) {
  // Release all profile counters
  uint64_t i;
  for (i=0;i<gem_profile.num_threads;++i) {
    mm_free(gem_profile.profile[i].timers);
    mm_free(gem_profile.profile[i].counters);
    mm_free(gem_profile.profile[i].ranks);
    stats_vector_delete(gem_profile.profile[i].strata_deltas_edit);
    stats_vector_delete(gem_profile.profile[i].strata_deltas_swg);
  }
  mm_free(gem_profile.profile);
}
/*
 * Accessors
 */
profile_t* PROF_GET_PROFILE() {
  return gem_profile.profile;
}
/*
 * PROFILE-TIME functions
 */
void PROF_START_TIMER(const uint64_t timer) {
  TIMER_START(gem_profile.profile[gtid()].timers+timer);
}
void PROF_STOP_TIMER(const uint64_t timer) {
  TIMER_STOP(gem_profile.profile[gtid()].timers+timer);
}
void PROF_PAUSE_TIMER(const uint64_t timer) {
  TIMER_PAUSE(gem_profile.profile[gtid()].timers+timer);
}
void PROF_CONTINUE_TIMER(const uint64_t timer) {
  TIMER_CONTINUE(gem_profile.profile[gtid()].timers+timer);
}
void PROF_RESET_TIMER(const uint64_t timer) {
  TIMER_RESET(gem_profile.profile[gtid()].timers+timer);
}
gem_timer_t* PROF_GET_TIMER(const uint64_t timer) {
  return gem_profile.profile[gtid()].timers+timer;
}
/*
 * PROFILE-COUNTERS functions
 */
void PROF_RESET_COUNTER(const uint64_t counter) {
  COUNTER_RESET(gem_profile.profile[gtid()].counters+counter);
}
void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value) {
  COUNTER_ADD(gem_profile.profile[gtid()].counters+counter,value);
}
void PROF_INC_COUNTER(const uint64_t counter) {
  COUNTER_ADD(gem_profile.profile[gtid()].counters+counter,1);
}
gem_counter_t* PROF_GET_COUNTER(const uint64_t counter) {
  return gem_profile.profile[gtid()].counters+counter;
}
/*
 * PROFILE-RANKS functions
 */
extern uint64_t _bwt_ranks; // Bwt rank counter
void PROF_START_RANK(const uint64_t rank) {
  RCOUNTER_START(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
void PROF_STOP_RANK(const uint64_t rank) {
  RCOUNTER_STOP(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
void PROF_PAUSE_RANK(const uint64_t rank) {
  RCOUNTER_PAUSE(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
void PROF_CONTINUE_RANK(const uint64_t rank) {
  RCOUNTER_CONTINUE(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
void PROF_RESET_RANK(const uint64_t rank) {
  RCOUNTER_RESET(gem_profile.profile[gtid()].ranks+rank);
}
gem_counter_t* PROF_GET_RANK(const uint64_t rank) {
  return &gem_profile.profile[gtid()].ranks[rank].counter;
}
/*
 * PROFILE-COMBINED (TIME/RANKS) functions
 */
void PROF_START(const uint64_t timer__ranks) {
  TIMER_START(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_START(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
void PROF_STOP(const uint64_t timer__ranks) {
  TIMER_STOP(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_STOP(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
void PROF_PAUSE(const uint64_t timer__ranks) {
  TIMER_PAUSE(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_PAUSE(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
void PROF_CONTINUE(const uint64_t timer__ranks) {
  TIMER_CONTINUE(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_CONTINUE(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
void PROF_RESET(const uint64_t timer__ranks) {
  TIMER_RESET(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_RESET(gem_profile.profile[gtid()].ranks+timer__ranks);
}
/*
 * Display statistics
 */
uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter) {
  const uint64_t gtid = gtid();
  return PERCENTAGE(
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters + counter),
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters + total_counter) );
}
double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2) {
  const uint64_t gtid = gtid();
  return DIV(
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters+counter1),
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters+counter2) );
}
double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t total_timer) {
  const uint64_t gtid = gtid();
  const uint64_t timer_ns = TIMER_GET_TOTAL_NS(gem_profile.profile[gtid].timers+timer);
  const uint64_t total_timer_ns = TIMER_GET_TOTAL_NS(gem_profile.profile[gtid].timers+total_timer);
  return PERCENTAGE(timer_ns,total_timer_ns);
}
double PROF_TIME_PER_CALL(const uint64_t timer) {
  const uint64_t gtid = gtid();
  gem_timer_t* const prof_timer = gem_profile.profile[gtid].timers+timer;
  const uint64_t num_calls = TIMER_GET_NUM_SAMPLES(prof_timer);
  const uint64_t timer_ns = TIMER_GET_TOTAL_NS(prof_timer);
  return (num_calls!=0) ? (TIMER_CONVERT_NS_TO_MS(timer_ns) / (double)num_calls) : 0.0;
}
/*
 * Utils
 */
void PROF_REDUCE_SUM(void) {
  uint64_t i;
  for (i=1;i<gem_profile.num_threads;++i) {
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      COUNTER_COMBINE_SUM(&gem_profile.profile[0].timers[j].time_ns,&gem_profile.profile[i].timers[j].time_ns);
      COUNTER_COMBINE_SUM(gem_profile.profile[0].counters+j,gem_profile.profile[i].counters+j);
      COUNTER_COMBINE_SUM(&gem_profile.profile[0].ranks[j].counter,&gem_profile.profile[i].ranks[j].counter);
    }
  }
}
void PROF_REDUCE_MAX(void) {
  uint64_t i;
  for (i=1;i<gem_profile.num_threads;++i) {
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      COUNTER_COMBINE_MAX(&gem_profile.profile[0].timers[j].time_ns,&gem_profile.profile[i].timers[j].time_ns);
      COUNTER_COMBINE_MAX(gem_profile.profile[0].counters+j,gem_profile.profile[i].counters+j);
      COUNTER_COMBINE_MAX(&gem_profile.profile[0].ranks[j].counter,&gem_profile.profile[i].ranks[j].counter);
    }
  }
}
void PROF_REDUCE_MIN(void) {
  uint64_t i;
  for (i=1;i<gem_profile.num_threads;++i) {
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      COUNTER_COMBINE_MIN(&gem_profile.profile[0].timers[j].time_ns,&gem_profile.profile[i].timers[j].time_ns);
      COUNTER_COMBINE_MIN(gem_profile.profile[0].counters+j,gem_profile.profile[i].counters+j);
      COUNTER_COMBINE_MIN(&gem_profile.profile[0].ranks[j].counter,&gem_profile.profile[i].ranks[j].counter);
    }
  }
}
void PROF_REDUCE_MEAN(void) {
  uint64_t i;
  for (i=1;i<gem_profile.num_threads;++i) {
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      COUNTER_COMBINE_MEAN(&gem_profile.profile[0].timers[j].time_ns,&gem_profile.profile[i].timers[j].time_ns);
      COUNTER_COMBINE_MEAN(gem_profile.profile[0].counters+j,gem_profile.profile[i].counters+j);
      COUNTER_COMBINE_MEAN(&gem_profile.profile[0].ranks[j].counter,&gem_profile.profile[i].ranks[j].counter);
    }
  }
}
void PROF_REDUCE_SAMPLE(void) {
  uint64_t i, j;
  for (j=0;j<GP_MAX_COUNTERS;++j) {
    i=0;
    while (i<gem_profile.num_threads-1 && gem_profile.profile[i].timers[j].time_ns.total==0) i++;
    gem_profile.profile[0].timers[j] = gem_profile.profile[i].timers[j];
    i=0;
    while (i<gem_profile.num_threads-1 && gem_profile.profile[i].counters[j].total==0) i++;
    gem_profile.profile[0].counters[j] = gem_profile.profile[i].counters[j];
    i=0;
    while (i<gem_profile.num_threads-1 && gem_profile.profile[i].ranks[j].accumulated==0) i++;
    gem_profile.profile[0].ranks[j] = gem_profile.profile[i].ranks[j];
  }
}

#endif /* GEM_PROFILE */

