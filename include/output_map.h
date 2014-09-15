/*
 * PROJECT: GEMMapper
 * FILE: output_map.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef OUTPUT_MAP_H_
#define OUTPUT_MAP_H_

#include "essentials.h"
#include "buffered_output_file.h"
#include "sequence.h"
#include "matches.h"

GEM_INLINE void output_map_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches);

#endif /* OUTPUT_MAP_H_ */
