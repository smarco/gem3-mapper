/*
 * PROJECT: GEMMapper
 * FILE: sequence_qualities_model.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Paolo Ribeca <paolo.ribeca@gmail.com>
 * DESCRIPTION:
 */

#ifndef SEQUENCE_QUALITIES_MODEL_H_
#define SEQUENCE_QUALITIES_MODEL_H_

#include "utils/essentials.h"
#include "text/sequence.h"

/*
 * Quality Data Structures
 */
typedef enum {
  sequence_qualities_model_flat,
  sequence_qualities_model_gem
} sequence_qualities_model_t;
typedef enum {
  sequence_qualities_ignore,
  sequence_qualities_offset_33,
  sequence_qualities_offset_64
} sequence_qualities_format_t;
typedef enum {
  qm_real=0,
  qm_pseudo=1
} quality_type_t;

/*
 * Quality Models
 */
void sequence_qualities_model_process(
    sequence_t* const sequence,
    const sequence_qualities_model_t qualities_model,
    const sequence_qualities_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask);

#endif /* SEQUENCE_QUALITIES_MODEL_H_ */
