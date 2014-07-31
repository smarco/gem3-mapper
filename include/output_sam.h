/*
 * PROJECT: GEMMapper
 * FILE: output_sam.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef OUTPUT_SAM_H_
#define OUTPUT_SAM_H_

#include "essentials.h"
#include "buffered_output_file.h"
#include "sequence.h"
#include "matches.h"

GEM_INLINE void output_sam_single_end_matches(
    output_buffer_t* const output_buffer,
    const sequence_t* const seq_read,const matches_t* const matches);

#endif /* OUTPUT_SAM_H_ */
