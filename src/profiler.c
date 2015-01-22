/*
 * PROJECT: GEMMapper
 * FILE: profiler.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple time/functional profiler module
 */

#include "profiler.h"
#include "gthread.h"

/*
 * Counters
 */
GEM_INLINE void COUNTER_RESET(gem_counter_t* const counter) {
  counter->total = 0.0;
  counter->samples = 0;
}
GEM_INLINE void COUNTER_ADD(gem_counter_t* const counter,const uint64_t amount) {
  // Add to total & increment number of samples
  counter->total += amount;
  ++(counter->samples);
  // See Knuth TAOCP vol 2, 3rd edition, page 232
  if (counter->samples == 1) {
    counter->min = amount;
    counter->max = amount;
    counter->m_oldM = amount;
    counter->m_newM = amount;
    counter->m_oldS = 0.0;
  } else {
    counter->min = MIN(counter->min,amount);
    counter->max = MAX(counter->max,amount);
    counter->m_newM = counter->m_oldM + ((double)amount-counter->m_oldM)/(double)counter->samples;
    counter->m_newS = counter->m_oldS + ((double)amount-counter->m_oldM)*((double)amount-counter->m_newM);
    counter->m_oldM = counter->m_newM;
    counter->m_oldS = counter->m_newS;
  }
}
GEM_INLINE uint64_t COUNTER_GET_TOTAL(const gem_counter_t* const counter) {
  return counter->total;
}
GEM_INLINE uint64_t COUNTER_GET_NUM_SAMPLES(const gem_counter_t* const counter) {
  return counter->samples;
}
GEM_INLINE uint64_t COUNTER_GET_MIN(const gem_counter_t* const counter) {
  return counter->min;
}
GEM_INLINE uint64_t COUNTER_GET_MAX(const gem_counter_t* const counter) {
  return counter->max;
}
GEM_INLINE double COUNTER_GET_MEAN(const gem_counter_t* const counter) {
  return (double)counter->total/(double)counter->samples;
}
GEM_INLINE double COUNTER_GET_VARIANCE(const gem_counter_t* const counter) {
  return ((counter->samples > 1) ? counter->m_newS/(double)(counter->samples - 1) : 0.0);
}
GEM_INLINE double COUNTER_GET_STDDEV(const gem_counter_t* const counter) {
  return sqrt(COUNTER_GET_VARIANCE(counter));
}
GEM_INLINE void COUNTER_COMBINE_SUM(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
  counter_dst->total += counter_src->total;
  counter_dst->samples += counter_src->samples;
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS; // FIXME
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
}
GEM_INLINE void COUNTER_COMBINE_MAX(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
  if (counter_dst->total < counter_src->total) {
    counter_dst->total = counter_src->total;
    counter_dst->samples = counter_src->samples;
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
  }
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
}
GEM_INLINE void COUNTER_COMBINE_MIN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
  if (counter_dst->samples==0 || counter_dst->total==0 || counter_dst->total > counter_src->total) {
    counter_dst->total = counter_src->total;
    counter_dst->samples = counter_src->samples;
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
  }
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
}
GEM_INLINE void COUNTER_COMBINE_MEAN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
  // FIXME Horrible, but listen! Now, I just want a roughly estimation
  counter_dst->total = (counter_dst->total+counter_src->total)/2;
  counter_dst->samples = (counter_dst->samples+counter_src->samples)/2;
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
}
GEM_INLINE void COUNTER_PRINT(
    FILE* const stream,const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,const char* const units,const bool full_report) {
  const uint64_t total = COUNTER_GET_TOTAL(counter);
  // Print Total
  if (total >= BUFFER_SIZE_1G) {
    fprintf(stream,"%7.2f G%s",(double)total/BUFFER_SIZE_1G,units);
  } else if (total >= BUFFER_SIZE_1M) {
    fprintf(stream,"%7.2f M%s ",(double)total/BUFFER_SIZE_1M,units);
  } else if (total >= BUFFER_SIZE_1K) {
    fprintf(stream,"%7.2f K%s",(double)total/BUFFER_SIZE_1K,units);
  } else {
    fprintf(stream,"%7.2f %s  ",(double)total,units);
  }
  // Print percentage wrt reference
  if (ref_counter!=NULL) {
    if (total==0) {
        fprintf(stream," (  0.00 %%)");
    } else {
      const uint64_t total_ref = COUNTER_GET_TOTAL(ref_counter);
      if (total_ref==0) {
        fprintf(stream," (  n/a  %%)");
      } else {
        const double percentage = (double)(total*100)/(double)total_ref;
        fprintf(stream," (%6.02f %%)",percentage);
      }
    }
  }
  // Full report
  if (!full_report) {
    fprintf(stream,"\n");
    return;
  }
  // Print Samples
  const uint64_t num_samples = COUNTER_GET_NUM_SAMPLES(counter);
  if (num_samples >= BUFFER_SIZE_1G) {
    fprintf(stream," (samples=%luG",num_samples/BUFFER_SIZE_1G);
  } else if (num_samples >= BUFFER_SIZE_1M) {
    fprintf(stream," (samples=%luM",num_samples/BUFFER_SIZE_1M);
  } else if (num_samples >= BUFFER_SIZE_1K) {
    fprintf(stream," (samples=%luK",num_samples/BUFFER_SIZE_1K);
  } else {
    fprintf(stream," (samples=%lu",num_samples);
    if (num_samples==0) {
      fprintf(stream,",--n/a--)}\n");
      return;
    }
  }
  // Print Mean
  const double mean = COUNTER_GET_MEAN(counter);
  if (mean >= BUFFER_SIZE_1G) {
    fprintf(stream,"{mean%.2fG",mean/BUFFER_SIZE_1G);
  } else if (mean >= BUFFER_SIZE_1M) {
    fprintf(stream,"{mean%.2fM",mean/BUFFER_SIZE_1M);
  } else if (mean >= BUFFER_SIZE_1K) {
    fprintf(stream,"{mean%.2fK",mean/BUFFER_SIZE_1K);
  } else {
    fprintf(stream,"{mean%.2f",mean);
  }
  // Print Min
  const uint64_t min = COUNTER_GET_MIN(counter);
  if (min >= BUFFER_SIZE_1G) {
    fprintf(stream,",min%.2fG",(double)min/BUFFER_SIZE_1G);
  } else if (min >= BUFFER_SIZE_1M) {
    fprintf(stream,",min%.2fM",(double)min/BUFFER_SIZE_1M);
  } else if (min >= BUFFER_SIZE_1K) {
    fprintf(stream,",min%.2fK",(double)min/BUFFER_SIZE_1K);
  } else {
    fprintf(stream,",min%.2f",(double)min);
  }
  // Print Max
  const uint64_t max = COUNTER_GET_MAX(counter);
  if (max >= BUFFER_SIZE_1G) {
    fprintf(stream,",Max%.2fG",(double)max/BUFFER_SIZE_1G);
  } else if (max >= BUFFER_SIZE_1M) {
    fprintf(stream,",Max%.2fM",(double)max/BUFFER_SIZE_1M);
  } else if (max >= BUFFER_SIZE_1K) {
    fprintf(stream,",Max%.2fK",(double)max/BUFFER_SIZE_1K);
  } else {
    fprintf(stream,",Max%.2f",(double)max);
  }
  // Print Variance
  const uint64_t var = COUNTER_GET_VARIANCE(counter);
  if (var >= BUFFER_SIZE_1G) {
    fprintf(stream,",Var%.2fG",(double)var/BUFFER_SIZE_1G);
  } else if (var >= BUFFER_SIZE_1M) {
    fprintf(stream,",Var%.2fM",(double)var/BUFFER_SIZE_1M);
  } else if (var >= BUFFER_SIZE_1K) {
    fprintf(stream,",Var%.2fK",(double)var/BUFFER_SIZE_1K);
  } else {
    fprintf(stream,",Var%.2f",(double)var);
  }
  // Print Standard Deviation
  const uint64_t stdDev = COUNTER_GET_STDDEV(counter);
  if (stdDev >= BUFFER_SIZE_1G) {
    fprintf(stream,",StdDev%.2fG)}\n",(double)stdDev/BUFFER_SIZE_1G);
  } else if (stdDev >= BUFFER_SIZE_1M) {
    fprintf(stream,",StdDev%.2fM)}\n",(double)stdDev/BUFFER_SIZE_1M);
  } else if (stdDev >= BUFFER_SIZE_1K) {
    fprintf(stream,",StdDev%.2fK)}\n",(double)stdDev/BUFFER_SIZE_1K);
  } else {
    fprintf(stream,",StdDev%.2f)}\n",(double)stdDev);
  }
}
GEM_INLINE void PERCENTAGE_PRINT(FILE* const stream,const gem_counter_t* const counter) {
  // Print Mean
  const double mean = COUNTER_GET_MEAN(counter);
  fprintf(stream,"%7.2f%%",mean);
  // Print Samples
  const uint64_t num_samples = COUNTER_GET_NUM_SAMPLES(counter);
  if (num_samples >= BUFFER_SIZE_1G) {
    fprintf(stream," (samples=%luG",num_samples/BUFFER_SIZE_1G);
  } else if (num_samples >= BUFFER_SIZE_1M) {
    fprintf(stream," (samples=%luM",num_samples/BUFFER_SIZE_1M);
  } else if (num_samples >= BUFFER_SIZE_1K) {
    fprintf(stream," (samples=%luK",num_samples/BUFFER_SIZE_1K);
  } else {
    fprintf(stream," (samples=%lu",num_samples);
  }
  if (num_samples == 0) {
    fprintf(stream,")\n");
    return;
  }
  // Print Min/Max
  fprintf(stream,",min%.2f%%,Max%.2f%%",(double)COUNTER_GET_MIN(counter),(double)COUNTER_GET_MAX(counter));
  // Print Variance/StandardDeviation
  fprintf(stream,",Var%.2f,StdDev%.2f)\n",COUNTER_GET_VARIANCE(counter),COUNTER_GET_STDDEV(counter));
}
/*
 * Timers
 */
GEM_INLINE void TIMER_START(gem_timer_t* const timer) {
  timer->accumulated = 0;
  TIMER_CONTINUE(timer);
}
GEM_INLINE void TIMER_STOP(gem_timer_t* const timer) {
  TIMER_PAUSE(timer);
  COUNTER_ADD(&timer->time_ns,timer->accumulated);
}
GEM_INLINE void TIMER_PAUSE(gem_timer_t* const timer) {
  clock_gettime(CLOCK_REALTIME,&timer->end_timer);
  timer->accumulated += TIME_DIFF_NS(timer->begin_timer,timer->end_timer);
}
GEM_INLINE void TIMER_CONTINUE(gem_timer_t* const timer) {
  clock_gettime(CLOCK_REALTIME,&timer->begin_timer);
}
GEM_INLINE void TIMER_RESET(gem_timer_t* const timer) {
  COUNTER_RESET(&timer->time_ns);
}
GEM_INLINE void TIMER_RESTART(gem_timer_t* const timer) {
  TIMER_RESET(timer);
  TIMER_START(timer);
}
GEM_INLINE uint64_t TIMER_GET_TOTAL_NS(const gem_timer_t* const timer) {
  return COUNTER_GET_TOTAL(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_NUM_SAMPLES(const gem_timer_t* const timer) {
  return COUNTER_GET_NUM_SAMPLES(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_MIN_NS(const gem_timer_t* const timer) {
  return COUNTER_GET_MIN(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_MAX_NS(const gem_timer_t* const timer) {
  return COUNTER_GET_MAX(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_MEAN(const gem_timer_t* const timer) {
  return COUNTER_GET_MEAN(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_VARIANCE(const gem_timer_t* const timer) {
  return COUNTER_GET_VARIANCE(&timer->time_ns);
}
GEM_INLINE uint64_t TIMER_GET_STDDEV(const gem_timer_t* const timer) {
  return COUNTER_GET_STDDEV(&timer->time_ns);
}
GEM_INLINE void TIMER_PRINT(
    FILE* const stream,const gem_timer_t* const timer,const gem_timer_t* const ref_timer) {
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
    fprintf(stream,"%7lu ns",total_time_ns);
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
    fprintf(stream," (%5lu Gcalls",num_calls/1000000000);
  } else if (num_calls >= 1000000) {
    fprintf(stream," (%5lu Mcalls",num_calls/1000000);
  } else if (num_calls >= 1000) {
    fprintf(stream," (%5lu Kcalls",num_calls/1000);
  } else if (num_calls > 1 || num_calls == 0) {
    fprintf(stream," (%5lu  calls",num_calls);
  } else {
    fprintf(stream," (%5lu   call",num_calls);
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
      fprintf(stream,",%7lu ns/call",ns_per_call);
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
    fprintf(stream," {min%luns",min_ns);
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
    fprintf(stream,",Max%luns})\n",max_ns);
  }
}
/*
 * Reference Counter (Counts wrt a reference counter. Eg ranks)
 */
GEM_INLINE void RCOUNTER_START(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  rcounter->accumulated = 0;
  RCOUNTER_CONTINUE(rcounter,reference);
}
GEM_INLINE void RCOUNTER_STOP(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  RCOUNTER_PAUSE(rcounter,reference);
  COUNTER_ADD(&rcounter->counter,rcounter->accumulated);
}
GEM_INLINE void RCOUNTER_PAUSE(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  rcounter->accumulated += reference-rcounter->begin_count;
}
GEM_INLINE void RCOUNTER_CONTINUE(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  rcounter->begin_count = reference;
}
GEM_INLINE void RCOUNTER_RESET(gem_reference_counter_t* const rcounter) {
  COUNTER_RESET(&rcounter->counter);
}
GEM_INLINE uint64_t RCOUNTER_GET_TOTAL(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_TOTAL(&rcounter->counter);
}
GEM_INLINE uint64_t RCOUNTER_GET_NUM_SAMPLES(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_NUM_SAMPLES(&rcounter->counter);
}
GEM_INLINE uint64_t RCOUNTER_GET_MIN(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_MIN(&rcounter->counter);
}
GEM_INLINE uint64_t RCOUNTER_GET_MAX(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_MAX(&rcounter->counter);
}
GEM_INLINE uint64_t RCOUNTER_GET_MEAN(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_MEAN(&rcounter->counter);
}
GEM_INLINE uint64_t RCOUNTER_GET_VARIANCE(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_VARIANCE(&rcounter->counter);
}
GEM_INLINE uint64_t RCOUNTER_GET_STDDEV(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_STDDEV(&rcounter->counter);
}
/*
 * Profile
 */
#ifndef GEM_NOPROFILE
typedef struct {
  gem_timer_t* timers;            // Time counters
  gem_counter_t* counters;        // General counters & functions
  gem_reference_counter_t* ranks; // Ranks counters
} profile_t;
typedef struct {
  // Profiler
  profile_t* profile;
  // Limits
  uint64_t num_threads;
} profiler_t;
profiler_t gem_profile; // THE GREAT PROFILER

/*
 * Setup
 */
GEM_INLINE void PROF_NEW(const uint64_t num_threads) {
  // Allocate handler
  gem_profile.profile = mm_calloc(num_threads,profile_t,true);
  // Initialize
  gem_profile.num_threads = num_threads;
  // Allocate all profiles
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    gem_profile.profile[i].timers = mm_calloc(GP_MAX_COUNTERS,gem_timer_t,true);
    gem_profile.profile[i].counters = mm_calloc(GP_MAX_COUNTERS,gem_counter_t,true);
    gem_profile.profile[i].ranks = mm_calloc(GP_MAX_COUNTERS,gem_reference_counter_t,true);
    // Initialize minimums (towards combine helps a lot)
    uint64_t j;
    for (j=0;j<GP_MAX_COUNTERS;++j) {
      gem_profile.profile[i].timers[j].time_ns.min = UINT64_MAX;
      gem_profile.profile[i].counters[j].min = UINT64_MAX;
      gem_profile.profile[i].ranks[j].counter.min = UINT64_MAX;
    }
  }
}
GEM_INLINE void PROF_DELETE() {
  // Release all profile counters
  uint64_t i;
  for (i=0;i<gem_profile.num_threads;++i) {
    mm_free(gem_profile.profile[i].timers);
    mm_free(gem_profile.profile[i].counters);
    mm_free(gem_profile.profile[i].ranks);
  }
  mm_free(gem_profile.profile);
}
/*
 * PROFILE-TIME functions
 */
GEM_INLINE void PROF_START_TIMER(const uint64_t timer) {
  TIMER_START(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_STOP_TIMER(const uint64_t timer) {
  TIMER_STOP(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_PAUSE_TIMER(const uint64_t timer) {
  TIMER_PAUSE(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_CONTINUE_TIMER(const uint64_t timer) {
  TIMER_CONTINUE(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE void PROF_RESET_TIMER(const uint64_t timer) {
  TIMER_RESET(gem_profile.profile[gtid()].timers+timer);
}
GEM_INLINE gem_timer_t* PROF_GET_TIMER(const uint64_t timer) {
  return gem_profile.profile[gtid()].timers+timer;
}
/*
 * PROFILE-COUNTERS functions
 */
GEM_INLINE void PROF_RESET_COUNTER(const uint64_t counter) {
  COUNTER_RESET(gem_profile.profile[gtid()].counters+counter);
}
GEM_INLINE void PROF_ADD_COUNTER(const uint64_t counter,const uint64_t value) {
  COUNTER_ADD(gem_profile.profile[gtid()].counters+counter,value);
}
GEM_INLINE void PROF_INC_COUNTER(const uint64_t counter) {
  COUNTER_ADD(gem_profile.profile[gtid()].counters+counter,1);
}
GEM_INLINE gem_counter_t* PROF_GET_COUNTER(const uint64_t counter) {
  return gem_profile.profile[gtid()].counters+counter;
}
/*
 * PROFILE-RANKS functions
 */
extern uint64_t _bwt_ranks;
GEM_INLINE void PROF_START_RANK(const uint64_t rank) {
  RCOUNTER_START(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
GEM_INLINE void PROF_STOP_RANK(const uint64_t rank) {
  RCOUNTER_STOP(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
GEM_INLINE void PROF_PAUSE_RANK(const uint64_t rank) {
  RCOUNTER_PAUSE(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
GEM_INLINE void PROF_CONTINUE_RANK(const uint64_t rank) {
  RCOUNTER_CONTINUE(gem_profile.profile[gtid()].ranks+rank,_bwt_ranks);
}
GEM_INLINE void PROF_RESET_RANK(const uint64_t rank) {
  RCOUNTER_RESET(gem_profile.profile[gtid()].ranks+rank);
}
GEM_INLINE gem_counter_t* PROF_GET_RANK(const uint64_t rank) {
  return &gem_profile.profile[gtid()].ranks[rank].counter;
}
/*
 * PROFILE-COMBINED (TIME/RANKS) functions
 */
GEM_INLINE void PROF_START(const uint64_t timer__ranks) {
  TIMER_START(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_START(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
GEM_INLINE void PROF_STOP(const uint64_t timer__ranks) {
  TIMER_STOP(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_STOP(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
GEM_INLINE void PROF_PAUSE(const uint64_t timer__ranks) {
  TIMER_PAUSE(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_PAUSE(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
GEM_INLINE void PROF_CONTINUE(const uint64_t timer__ranks) {
  TIMER_CONTINUE(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_CONTINUE(gem_profile.profile[gtid()].ranks+timer__ranks,_bwt_ranks);
}
GEM_INLINE void PROF_RESET(const uint64_t timer__ranks) {
  TIMER_RESET(gem_profile.profile[gtid()].timers+timer__ranks);
  RCOUNTER_RESET(gem_profile.profile[gtid()].ranks+timer__ranks);
}
/*
 * Display statistics
 */
GEM_INLINE uint64_t PROF_COUNT_PERCENTAGE(const uint64_t counter,const uint64_t total_counter) {
  const uint64_t gtid = gtid();
  return PERCENTAGE(
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters + counter),
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters + total_counter) );
}
GEM_INLINE double PROF_COUNT_DIV(const uint64_t counter1,const uint64_t counter2) {
  const uint64_t gtid = gtid();
  return DIV(
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters+counter1),
      COUNTER_GET_TOTAL(gem_profile.profile[gtid].counters+counter2) );
}
GEM_INLINE double PROF_TIME_PERCENTAGE(const uint64_t timer,const uint64_t total_timer) {
  const uint64_t gtid = gtid();
  const uint64_t timer_ns = TIMER_GET_TOTAL_NS(gem_profile.profile[gtid].timers+timer);
  const uint64_t total_timer_ns = TIMER_GET_TOTAL_NS(gem_profile.profile[gtid].timers+total_timer);
  return PERCENTAGE(timer_ns,total_timer_ns);
}
GEM_INLINE double PROF_TIME_PER_CALL(const uint64_t timer) {
  const uint64_t gtid = gtid();
  gem_timer_t* const prof_timer = gem_profile.profile[gtid].timers+timer;
  const uint64_t num_calls = TIMER_GET_NUM_SAMPLES(prof_timer);
  const uint64_t timer_ns = TIMER_GET_TOTAL_NS(prof_timer);
  return (num_calls!=0) ? (TIMER_CONVERT_NS_TO_MS(timer_ns) / (double)num_calls) : 0.0;
}
/*
 * Utils
 */
GEM_INLINE void PROF_REDUCE_SUM() {
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
GEM_INLINE void PROF_REDUCE_MAX() {
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
GEM_INLINE void PROF_REDUCE_MIN() {
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
GEM_INLINE void PROF_REDUCE_MEAN() {
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
GEM_INLINE void PROF_REDUCE_SAMPLE() {
  uint64_t i, j;
  for (j=0;j<GP_MAX_COUNTERS;++j) {
    i=0;
    while (i<gem_profile.num_threads-1 && gem_profile.profile[i].timers[j].time_ns.total==0) i++;
    gem_profile.profile[0].timers[j].time_ns = gem_profile.profile[i].timers[j].time_ns;
    i=0;
    while (i<gem_profile.num_threads-1 && gem_profile.profile[i].counters[j].total==0) i++;
    gem_profile.profile[0].counters[j] = gem_profile.profile[i].counters[j];
    i=0;
    while (i<gem_profile.num_threads-1 && gem_profile.profile[i].ranks[j].accumulated==0) i++;
    gem_profile.profile[0].ranks[j] = gem_profile.profile[i].ranks[j];
  }
}
#endif /* !GEM_NOPROFILE */

