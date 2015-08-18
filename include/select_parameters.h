/*
 * PROJECT: GEMMapper
 * FILE: select_parameters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SELECT_PARAMETERS_H_
#define SELECT_PARAMETERS_H_

#include "essentials.h"

/*
 * Select Parameters
 */
typedef enum {
  mapq_model_none,     // None
  mapq_model_gem,      // GEM Score (Case stratification + Logistic Regression)
  mapq_model_classify, // GEM Classification towards score calibration
} mapq_model_t;
typedef enum {
  matches_sorting_distance,
  matches_sorting_mapq
} matches_sorting_t;
typedef struct {
  /* MAPQ Score */
  mapq_model_t mapq_model;
  uint8_t mapq_threshold;
  /* Reporting */
  double min_reported_strata;
  uint64_t min_reported_strata_nominal;
  uint64_t min_reported_matches;
  uint64_t max_reported_matches;
  /* Check */
  bool check_correct;
  bool check_optimum;
  bool check_complete;
} select_parameters_t;

/*
 * Select Parameters Setup
 */
GEM_INLINE void select_parameters_init(select_parameters_t* const select_parameters);

GEM_INLINE void select_configure_reporting(
    select_parameters_t* const select_parameters,const float min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches);

GEM_INLINE void select_instantiate_values(
    select_parameters_t* const select_parameters,const uint64_t sequence_length);

#endif /* SELECT_PARAMETERS_H_ */
