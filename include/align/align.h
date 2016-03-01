/*
 * PROJECT: GEMMapper
 * FILE: align.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ALIGN_H_
#define ALIGN_H_

#include "utils/essentials.h"

/*
 * Constants
 */
#define ALIGN_DISTANCE_INF     UINT32_MAX
#define ALIGN_COLUMN_INF       UINT64_MAX
#define ALIGN_DISABLED         (UINT32_MAX-1)

/*
 * Check matches (CIGAR string against text & pattern)
 */
bool align_check(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    uint64_t const cigar_offset,
    uint64_t const cigar_length,
    const bool verbose);

/*
 * Compute edit distance (Basic DP-Matrix Alignment)
 */
int64_t align_dp_compute_edit_distance(
    const char* const key,
    const uint64_t key_length,
    const char* const text,
    const uint64_t text_length,
    const bool ends_free,
    uint64_t* const position);

#endif /* ALIGN_H_ */
