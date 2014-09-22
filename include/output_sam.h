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
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches);

#endif /* OUTPUT_SAM_H_ */
