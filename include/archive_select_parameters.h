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
#define check_none        0x0ull  /* Check nothing */
#define check_correctness 0x1ull  /* Check that the reported mappings are correct (position+CIGAR)*/
#define check_optimum     0x2ull  /* Check that the reported CIGAR optimizes the aligment-model */
#define check_completness 0x4ull  /* Check that within the MCS reported, no match is missing */
typedef enum {
  mapq_model_none,
  mapq_model_li,
  mapq_model_heath,
  mapq_model_sm,
  mapq_model_pr
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
  uint64_t check_matches_mask;
} select_parameters_t;

/*
 * Select Parameters Setup
 */
GEM_INLINE void archive_select_parameters_init(select_parameters_t* const select_parameters);

GEM_INLINE void archive_select_configure_reporting(
    select_parameters_t* const select_parameters,
    float min_decoded_strata,uint64_t max_decoded_matches,
    uint64_t min_reported_matches,uint64_t max_reported_matches);
GEM_INLINE void archive_select_configure_check(
    select_parameters_t* const select_parameters,const uint64_t check_matches_mask);

GEM_INLINE void archive_select_instantiate_values(
    select_parameters_t* const select_parameters,const uint64_t sequence_length);
/*
 * Error Messages
 */
//#define GEM_ERROR_ARCHIVE_SELECT_PARAMETERS_

#endif /* ARCHIVE_SELECT_PARAMETERS_H_ */
