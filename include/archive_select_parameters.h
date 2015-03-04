/*
 * PROJECT: GEMMapper
 * FILE: archive_select_parameters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SELECT_PARAMETERS_H_
#define ARCHIVE_SELECT_PARAMETERS_H_

#include "essentials.h"

/*
 * Select Parameters
 */
typedef enum {
  mapq_model_none,
  mapq_model_exp_relative_distance,
  mapq_model_exp_relative_score,
  mapq_model_test
} mapq_model_t;
typedef struct {
  /* MAPQ Score */
  mapq_model_t mapq_model;
  /* Reporting */
  double min_decoded_strata;
  uint64_t min_decoded_strata_nominal;
  uint64_t max_decoded_matches;
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
GEM_INLINE void archive_select_parameters_init(select_parameters_t* const select_parameters);

GEM_INLINE void archive_select_configure_reporting(
    select_parameters_t* const select_parameters,
    float min_decoded_strata,uint64_t max_decoded_matches,
    uint64_t min_reported_matches,uint64_t max_reported_matches);

GEM_INLINE void archive_select_instantiate_values(
    select_parameters_t* const select_parameters,const uint64_t sequence_length);
/*
 * Error Messages
 */
//#define GEM_ERROR_ARCHIVE_SELECT_PARAMETERS_

#endif /* ARCHIVE_SELECT_PARAMETERS_H_ */
