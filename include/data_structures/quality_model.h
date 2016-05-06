/*
 * PROJECT: GEMMapper
 * FILE: quality_model.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Paolo Ribeca <paolo.ribeca@gmail.com>
 * DESCRIPTION:
 */

#ifndef QUALITY_MODEL_H_
#define QUALITY_MODEL_H_

#include "utils/essentials.h"
#include "data_structures/sequence.h"

/*
 * Quality Data Structures
 */
typedef enum { quality_model_type_flat, quality_model_type_gem } quality_model_t;
typedef enum { qualities_ignore, qualities_offset_33, qualities_offset_64 } quality_format_t;

typedef enum { qm_real=0, qm_pseudo=1 } quality_type_t;

/*
 * Quality Model
 */
void quality_model(
    sequence_t* const sequence,
    const quality_model_t quality_model,
    const quality_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask);
void quality_model_flat(
    sequence_t* const sequence,
    const quality_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask);
void quality_model_gem(
    sequence_t* const sequence,
    const quality_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask);

#endif /* QUALITY_MODEL_H_ */
