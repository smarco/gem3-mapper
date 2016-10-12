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

#include "profiler/profiler_timer.h"

/*
 * Timers
 */
void TIMER_START(gem_timer_t* const timer) {
  timer->accumulated = 0;
  TIMER_CONTINUE(timer);
}
void TIMER_STOP(gem_timer_t* const timer) {
  TIMER_PAUSE(timer);
  COUNTER_ADD(&timer->time_ns,timer->accumulated);
}
void TIMER_PAUSE(gem_timer_t* const timer) {
  system_get_time(&timer->end_timer);
  timer->accumulated += TIME_DIFF_NS(timer->begin_timer,timer->end_timer);
}
void TIMER_CONTINUE(gem_timer_t* const timer) {
  system_get_time(&timer->begin_timer);
}
void TIMER_RESET(gem_timer_t* const timer) {
  COUNTER_RESET(&timer->time_ns);
}
void TIMER_RESTART(gem_timer_t* const timer) {
  TIMER_RESET(timer);
  TIMER_START(timer);
}
uint64_t TIMER_GET_TOTAL_NS(const gem_timer_t* const timer) {
  return COUNTER_GET_TOTAL(&timer->time_ns);
}
uint64_t TIMER_GET_NUM_SAMPLES(const gem_timer_t* const timer) {
  return COUNTER_GET_NUM_SAMPLES(&timer->time_ns);
}
uint64_t TIMER_GET_MIN_NS(const gem_timer_t* const timer) {
  return COUNTER_GET_MIN(&timer->time_ns);
}
uint64_t TIMER_GET_MAX_NS(const gem_timer_t* const timer) {
  return COUNTER_GET_MAX(&timer->time_ns);
}
uint64_t TIMER_GET_MEAN(const gem_timer_t* const timer) {
  return COUNTER_GET_MEAN(&timer->time_ns);
}
uint64_t TIMER_GET_VARIANCE(const gem_timer_t* const timer) {
  return COUNTER_GET_VARIANCE(&timer->time_ns);
}
uint64_t TIMER_GET_STDDEV(const gem_timer_t* const timer) {
  return COUNTER_GET_STDDEV(&timer->time_ns);
}
void TIMER_PRINT(
    FILE* const stream,const gem_timer_t* const timer,
    const gem_timer_t* const ref_timer) {
  const uint64_t total_time_ns = TIMER_GET_TOTAL_NS(timer);
  // Print Total
  if (total_time_ns >= 60000000000ull) {
    fprintf(stream,"%7.2f m ",TIMER_CONVERT_NS_TO_M(total_time_ns));
  } else if (total_time_ns >= 1000000000) {
    fprintf(stream,"%7.2f s ",TIMER_CONVERT_NS_TO_S(total_time_ns));
  } else if (total_time_ns >= 1000000) {
    fprintf(stream,"%7.2f ms",TIMER_CONVERT_NS_TO_MS(total_time_ns));
  } else if (total_time_ns >= 1000) {
    fprintf(stream,"%7.2f us",TIMER_CONVERT_NS_TO_US(total_time_ns));
  } else {
    fprintf(stream,"%7"PRIu64" ns",total_time_ns);
  }
  // Print percentage wrt reference
  if (ref_timer!=NULL) {
    if (total_time_ns==0) {
        fprintf(stream," (  0.00 %%)");
    } else {
      const uint64_t total_ref_time_ns = TIMER_GET_TOTAL_NS(ref_timer);
      if (total_ref_time_ns==0) {
        fprintf(stream," (  n/a  %%)");
      } else {
        const double percentage = (double)(total_time_ns*100)/(double)total_ref_time_ns;
        fprintf(stream," (%6.02f %%)",percentage);
      }
    }
  }
  // Print Calls
  const uint64_t num_calls = TIMER_GET_NUM_SAMPLES(timer);
  if (num_calls >= 1000000000) {
    fprintf(stream," (%5"PRIu64" Gcalls",num_calls/1000000000);
  } else if (num_calls >= 1000000) {
    fprintf(stream," (%5"PRIu64" Mcalls",num_calls/1000000);
  } else if (num_calls >= 1000) {
    fprintf(stream," (%5"PRIu64" Kcalls",num_calls/1000);
  } else if (num_calls > 1 || num_calls == 0) {
    fprintf(stream," (%5"PRIu64"  calls",num_calls);
  } else {
    fprintf(stream," (%5"PRIu64"   call",num_calls);
  }
  // Print time/call
  if (num_calls==0) {
    fprintf(stream,",   n/a   s/call)\n");
    return;
  } else {
    const uint64_t ns_per_call = total_time_ns / num_calls;
    if (ns_per_call > 1000000000) {
      fprintf(stream,",%7.2f  s/call",TIMER_CONVERT_NS_TO_S(ns_per_call));
    } else if (ns_per_call > 1000000) {
      fprintf(stream,",%7.2f ms/call",TIMER_CONVERT_NS_TO_MS(ns_per_call));
    } else if (ns_per_call > 1000) {
      fprintf(stream,",%7.2f us/call",TIMER_CONVERT_NS_TO_US(ns_per_call));
    } else {
      fprintf(stream,",%7"PRIu64" ns/call",ns_per_call);
    }
  }
  // Print Max
  const uint64_t min_ns = TIMER_GET_MIN_NS(timer);
  if (min_ns > 1000000000) {
    fprintf(stream," {min%.2fs",TIMER_CONVERT_NS_TO_S(min_ns));
  } else if (min_ns > 1000000) {
    fprintf(stream," {min%.2fms",TIMER_CONVERT_NS_TO_MS(min_ns));
  } else if (min_ns > 1000) {
    fprintf(stream," {min%.2fus",TIMER_CONVERT_NS_TO_US(min_ns));
  } else {
    fprintf(stream," {min%"PRIu64"ns",min_ns);
  }
  // Print Min
  const uint64_t max_ns = TIMER_GET_MAX_NS(timer);
  if (max_ns > 1000000000) {
    fprintf(stream,",Max%.2fs})\n",TIMER_CONVERT_NS_TO_S(max_ns));
  } else if (max_ns > 1000000) {
    fprintf(stream,",Max%.2fms})\n",TIMER_CONVERT_NS_TO_MS(max_ns));
  } else if (max_ns > 1000) {
    fprintf(stream,",Max%.2fus})\n",TIMER_CONVERT_NS_TO_US(max_ns));
  } else {
    fprintf(stream,",Max%"PRIu64"ns})\n",max_ns);
  }
}
