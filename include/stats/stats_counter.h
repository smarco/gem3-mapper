/*
 * PROJECT: GEMMapper
 * FILE: stats_counter.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#ifndef STATS_COUNTER_H_
#define STATS_COUNTER_H_

#include "system/commons.h"

/*
 * Counters
 *   (from http://www.johndcook.com/standard_deviation.html)
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

void COUNTER_COMBINE_SUM(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);
void COUNTER_COMBINE_MAX(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);
void COUNTER_COMBINE_MIN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);
void COUNTER_COMBINE_MEAN(gem_counter_t* const counter_dst,gem_counter_t* const counter_src);

void COUNTER_PRINT(
    FILE* const stream,
    const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,
    const char* const units,
    const bool full_report);
void SAMPLER_PRINT(
    FILE* const stream,
    const gem_counter_t* const counter,
    const gem_counter_t* const ref_counter,
    const char* const units);
void PERCENTAGE_PRINT(
    FILE* const stream,
    const gem_counter_t* const counter,
    const char* const units);

/*
 * Reference Counter (Counts wrt a reference counter. Eg ranks)
 */
typedef struct {
  uint64_t begin_count;  // Counter
  gem_counter_t counter; // Total count & samples taken
  uint64_t accumulated;  // Total accumulated
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
 * Display
 */
#define PRIcounter "lu(#%"PRIu64",m%"PRIu64",M%"PRIu64",{%.2f})"
#define PRIcounterVal(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter),COUNTER_GET_MEAN(counter)
#define PRIcounterX "lu(#%"PRIu64",m%"PRIu64",M%"PRIu64",{%.2f,%.2f,%.2f})"
#define PRIcounterXVal(counter) \
  COUNTER_GET_TOTAL(counter),COUNTER_GET_NUM_SAMPLES(counter), \
  COUNTER_GET_MIN(counter),COUNTER_GET_MAX(counter), \
  COUNTER_GET_MEAN(counter),COUNTER_GET_VARIANCE(counter),COUNTER_GET_STDDEV(counter)

#endif /* STATS_COUNTER_H_ */
