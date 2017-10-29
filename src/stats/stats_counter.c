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

#include "stats/stats_counter.h"
#include "system/mm.h"

/*
 * Counters
 */
void COUNTER_RESET(gem_counter_t* const counter) {
  memset(counter,0,sizeof(gem_counter_t));
}
void COUNTER_ADD(gem_counter_t* const counter,const uint64_t amount) {
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
uint64_t COUNTER_GET_TOTAL(const gem_counter_t* const counter) {
  return counter->total;
}
uint64_t COUNTER_GET_NUM_SAMPLES(const gem_counter_t* const counter) {
  return counter->samples;
}
uint64_t COUNTER_GET_MIN(const gem_counter_t* const counter) {
  return counter->min;
}
uint64_t COUNTER_GET_MAX(const gem_counter_t* const counter) {
  return counter->max;
}
double COUNTER_GET_MEAN(const gem_counter_t* const counter) {
  return (double)counter->total/(double)counter->samples;
}
double COUNTER_GET_VARIANCE(const gem_counter_t* const counter) {
  return ((counter->samples > 1) ? counter->m_newS/(double)(counter->samples - 1) : 0.0);
}
double COUNTER_GET_STDDEV(const gem_counter_t* const counter) {
  return sqrt(COUNTER_GET_VARIANCE(counter));
}
void COUNTER_COMBINE_SUM(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
  counter_dst->total += counter_src->total;
  counter_dst->samples += counter_src->samples;
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
}
void COUNTER_COMBINE_MAX(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
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
void COUNTER_COMBINE_MIN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
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
void COUNTER_COMBINE_MEAN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src) {
  // FIXME Horrible, but listen! Now, I just want a rough estimation
  counter_dst->total = (counter_dst->total+counter_src->total)/2;
  counter_dst->samples = (counter_dst->samples+counter_src->samples)/2;
  counter_dst->min = MIN(counter_dst->min,counter_src->min);
  counter_dst->max = MAX(counter_dst->max,counter_src->max);
  if (counter_src->m_newS!=0.0) counter_dst->m_newS = counter_src->m_newS;
  if (counter_src->m_newM!=0.0) counter_dst->m_newM = counter_src->m_newM;
  if (counter_src->m_oldS!=0.0) counter_dst->m_oldS = counter_src->m_oldS;
  if (counter_src->m_oldM!=0.0) counter_dst->m_oldM = counter_src->m_oldM;
}
void COUNTER_PRINT_STATS(
    FILE* const stream,
    const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,
    const char* const units) {
  // Print Samples
  const uint64_t num_samples = COUNTER_GET_NUM_SAMPLES(counter);
  if (num_samples >= METRIC_FACTOR_1G) {
    fprintf(stream," (samples=%"PRIu64"G",num_samples/METRIC_FACTOR_1G);
  } else if (num_samples >= METRIC_FACTOR_1M) {
    fprintf(stream," (samples=%"PRIu64"M",num_samples/METRIC_FACTOR_1M);
  } else if (num_samples >= METRIC_FACTOR_1K) {
    fprintf(stream," (samples=%"PRIu64"K",num_samples/METRIC_FACTOR_1K);
  } else {
    fprintf(stream," (samples=%"PRIu64"",num_samples);
    if (num_samples==0) {
      fprintf(stream,",--n/a--)}\n");
      return;
    }
  }
  // Print Mean
  const double mean = COUNTER_GET_MEAN(counter);
  if (mean >= METRIC_FACTOR_1G) {
    fprintf(stream,"{mean%.2fG",mean/METRIC_FACTOR_1G);
  } else if (mean >= METRIC_FACTOR_1M) {
    fprintf(stream,"{mean%.2fM",mean/METRIC_FACTOR_1M);
  } else if (mean >= METRIC_FACTOR_1K) {
    fprintf(stream,"{mean%.2fK",mean/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,"{mean%.2f",mean);
  }
  // Print Min
  const uint64_t min = COUNTER_GET_MIN(counter);
  if (min >= METRIC_FACTOR_1G) {
    fprintf(stream,",min%.2fG",(double)min/METRIC_FACTOR_1G);
  } else if (min >= METRIC_FACTOR_1M) {
    fprintf(stream,",min%.2fM",(double)min/METRIC_FACTOR_1M);
  } else if (min >= METRIC_FACTOR_1K) {
    fprintf(stream,",min%.2fK",(double)min/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",min%.2f",(double)min);
  }
  // Print Max
  const uint64_t max = COUNTER_GET_MAX(counter);
  if (max >= METRIC_FACTOR_1G) {
    fprintf(stream,",Max%.2fG",(double)max/METRIC_FACTOR_1G);
  } else if (max >= METRIC_FACTOR_1M) {
    fprintf(stream,",Max%.2fM",(double)max/METRIC_FACTOR_1M);
  } else if (max >= METRIC_FACTOR_1K) {
    fprintf(stream,",Max%.2fK",(double)max/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",Max%.2f",(double)max);
  }
  // Print Variance
  const uint64_t var = COUNTER_GET_VARIANCE(counter);
  if (var >= METRIC_FACTOR_1G) {
    fprintf(stream,",Var%.2fG",(double)var/METRIC_FACTOR_1G);
  } else if (var >= METRIC_FACTOR_1M) {
    fprintf(stream,",Var%.2fM",(double)var/METRIC_FACTOR_1M);
  } else if (var >= METRIC_FACTOR_1K) {
    fprintf(stream,",Var%.2fK",(double)var/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",Var%.2f",(double)var);
  }
  // Print Standard Deviation
  const uint64_t stdDev = COUNTER_GET_STDDEV(counter);
  if (stdDev >= METRIC_FACTOR_1G) {
    fprintf(stream,",StdDev%.2fG)}\n",(double)stdDev/METRIC_FACTOR_1G);
  } else if (stdDev >= METRIC_FACTOR_1M) {
    fprintf(stream,",StdDev%.2fM)}\n",(double)stdDev/METRIC_FACTOR_1M);
  } else if (stdDev >= METRIC_FACTOR_1K) {
    fprintf(stream,",StdDev%.2fK)}\n",(double)stdDev/METRIC_FACTOR_1K);
  } else {
    fprintf(stream,",StdDev%.2f)}\n",(double)stdDev);
  }
}
void COUNTER_PRINT(
    FILE* const stream,
    const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,
    const char* const units,
    const bool full_report) {
  const uint64_t total = COUNTER_GET_TOTAL(counter);
  // Print Total
  if (total >= METRIC_FACTOR_1G) {
    fprintf(stream,"%7.2f G%s",(double)total/METRIC_FACTOR_1G,units);
  } else if (total >= METRIC_FACTOR_1M) {
    fprintf(stream,"%7.2f M%s",(double)total/METRIC_FACTOR_1M,units);
  } else if (total >= METRIC_FACTOR_1K) {
    fprintf(stream,"%7.2f K%s",(double)total/METRIC_FACTOR_1K,units);
  } else {
    fprintf(stream,"%7.2f %s ",(double)total,units);
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
  } else {
    fprintf(stream,"           ");
  }
  // Full report
  if (!full_report) {
    fprintf(stream,"\n");
    return;
  } else {
    COUNTER_PRINT_STATS(stream,counter,ref_counter,units);
  }
}
void SAMPLER_PRINT(
    FILE* const stream,
    const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,
    const char* const units) {
  fprintf(stream,"\t\t\t\t");
  COUNTER_PRINT_STATS(stream,counter,ref_counter,units);
}
void PERCENTAGE_PRINT(
    FILE* const stream,
    const gem_counter_t* const counter,
    const char* const units) {
  // Print Mean
  const double mean = COUNTER_GET_MEAN(counter);
  fprintf(stream,"%7.2f %%%s\t\t",mean,units);
  // Print Samples
  const uint64_t num_samples = COUNTER_GET_NUM_SAMPLES(counter);
  if (num_samples >= METRIC_FACTOR_1G) {
    fprintf(stream," (samples=%"PRIu64"G",num_samples/METRIC_FACTOR_1G);
  } else if (num_samples >= METRIC_FACTOR_1M) {
    fprintf(stream," (samples=%"PRIu64"M",num_samples/METRIC_FACTOR_1M);
  } else if (num_samples >= METRIC_FACTOR_1K) {
    fprintf(stream," (samples=%"PRIu64"K",num_samples/METRIC_FACTOR_1K);
  } else {
    fprintf(stream," (samples=%"PRIu64"",num_samples);
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
 * Reference Counter (Counts wrt a reference counter. Eg ranks)
 */
void RCOUNTER_START(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  rcounter->accumulated = 0;
  RCOUNTER_CONTINUE(rcounter,reference);
}
void RCOUNTER_STOP(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  RCOUNTER_PAUSE(rcounter,reference);
  COUNTER_ADD(&rcounter->counter,rcounter->accumulated);
}
void RCOUNTER_PAUSE(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  rcounter->accumulated += reference-rcounter->begin_count;
}
void RCOUNTER_CONTINUE(gem_reference_counter_t* const rcounter,const uint64_t reference) {
  rcounter->begin_count = reference;
}
void RCOUNTER_RESET(gem_reference_counter_t* const rcounter) {
  COUNTER_RESET(&rcounter->counter);
}
uint64_t RCOUNTER_GET_TOTAL(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_TOTAL(&rcounter->counter);
}
uint64_t RCOUNTER_GET_NUM_SAMPLES(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_NUM_SAMPLES(&rcounter->counter);
}
uint64_t RCOUNTER_GET_MIN(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_MIN(&rcounter->counter);
}
uint64_t RCOUNTER_GET_MAX(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_MAX(&rcounter->counter);
}
uint64_t RCOUNTER_GET_MEAN(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_MEAN(&rcounter->counter);
}
uint64_t RCOUNTER_GET_VARIANCE(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_VARIANCE(&rcounter->counter);
}
uint64_t RCOUNTER_GET_STDDEV(gem_reference_counter_t* const rcounter) {
  return COUNTER_GET_STDDEV(&rcounter->counter);
}
