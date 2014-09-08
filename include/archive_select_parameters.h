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
  alignment_model_none,
  alignment_model_hamming,
  alignment_model_levenshtein,
  alignment_model_gap_affine
} alignment_model_t;
typedef struct {
  /* Single-end Alignment */
  /* Paired-end Alignment */
  /* Alignment Score */
  alignment_model_t alignment_model;
  //  float model_hamming_max_distance; // FIXME -> Maybe curation removes all this before existing
  //  float model_levenshtein_max_distance;
  //  float model_gap_affine_max_distance;
  uint64_t matching_score;
  uint64_t mismatch_penalty;
  uint64_t gap_open_penalty;
  uint64_t gap_extension_penalty;
  /* Mapping Quality */
  /* Reporting */
  float min_decoded_strata;
  uint64_t min_decoded_strata_nominal;
  uint64_t max_decoded_matches;
  uint64_t min_reported_matches;
  uint64_t max_reported_matches;
} select_parameters_t;

/*
 * Select Parameters Setup
 */
GEM_INLINE void archive_select_parameters_init(select_parameters_t* const select_parameters);

GEM_INLINE void archive_select_configure_reporting(
    select_parameters_t* const select_parameters,
    float min_decoded_strata,uint64_t max_decoded_matches,
    uint64_t min_reported_matches,uint64_t max_reported_matches);

GEM_INLINE void archive_select_configure_alignment_model(
    select_parameters_t* const select_parameters,const alignment_model_t alignment_model);

GEM_INLINE void archive_select_instantiate_values(
    select_parameters_t* const select_parameters,const uint64_t sequence_length);
/*
 * Error Messages
 */
//#define GEM_ERROR_ARCHIVE_SELECT_PARAMETERS_

#endif /* ARCHIVE_SELECT_PARAMETERS_H_ */
