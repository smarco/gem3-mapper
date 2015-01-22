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
#include "paired_matches.h"

/*
 * Utils
 */
GEM_INLINE void output_map_cigar(
    FILE* const stream,match_trace_t* const match_trace,matches_t* const matches);
GEM_INLINE void output_map_alignment_pretty(
    FILE* const stream,match_trace_t* const match_trace,matches_t* const matches,
    uint8_t* const key,const uint64_t key_length,uint8_t* const text,
    const uint64_t text_length,mm_stack_t* const mm_stack);

/*
 * MAP SingleEnd
 */
GEM_INLINE void output_map_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches);

/*
 * MAP PairedEnd
 */
GEM_INLINE void output_map_paired_end_matches(
    buffered_output_file_t* const buffered_output_file,sequence_t* const seq_read_end1,
    sequence_t* const seq_read_end2,paired_matches_t* const paired_matches);

#endif /* OUTPUT_MAP_H_ */
