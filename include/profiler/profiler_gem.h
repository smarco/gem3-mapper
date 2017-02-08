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

#ifndef PROFILER_GEM_H_
#define PROFILER_GEM_H_

#include "system/commons.h"
#include "system/mm.h"
#include "profiler/profiler_timer.h"
#include "stats/stats_counter.h"
#include "stats/stats_vector.h"

/*
 * Profile
 */
typedef enum {
  reduce_sum,
  reduce_max,
  reduce_min,
  reduce_mean,
  reduce_sample
} profile_reduce_type;
typedef struct {
  /* General timers,counters & ranks */
  gem_timer_t* timers;                  // Time counters
  gem_counter_t* counters;              // General counters & functions
  gem_reference_counter_t* ranks;       // Ranks counters
  /* Other counters */
  stats_vector_t* strata_deltas_edit;   // Strata-deltas stats
  stats_vector_t* strata_deltas_swg;    // Strata-deltas stats
} profile_t;
typedef struct {
  // Profiler
  profile_t* profile;
  // Limits
  uint64_t num_threads;
} profiler_t;

#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */

/*
 * Profiler Types
 */
#define GP_MAX_COUNTERS 1000

/*
 * Profile Block
 */
#define PROF_BLOCK()

/*
 * Setup
 */
void PROF_NEW(const uint64_t num_threads);
void PROF_DELETE(void);

/*
 * Accessors
 */
profile_t* PROF_GET_PROFILE();

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
void PROF_START(const uint64_t timer__ranks);
void PROF_STOP(const uint64_t timer__ranks);
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
void PROF_REDUCE_SUM(void);
void PROF_REDUCE_MAX(void);
void PROF_REDUCE_MIN(void);
void PROF_REDUCE_MEAN(void);
void PROF_REDUCE_SAMPLE(void);

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
#endif /* GEM_PROFILE */
#endif /* PROFILER_GEM_H_ */
