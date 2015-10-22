/*
 * PROJECT: GEMMapper
 * FILE: matches_cigar.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure to store alignment matches {sequence,position,strand,...}
 */

#ifndef MATCHES_CIGAR_H_
#define MATCHES_CIGAR_H_

#include "essentials.h"
#include "locator.h"
#include "interval_set.h"
#include "text_collection.h"
#include "match_elements.h"
#include "matches_counters.h"
#include "matches_metrics.h"
#include "align_swg.h"

/*
 * CIGAR Buffer Handling
 */
void matches_cigar_buffer_add_cigar_element(
    cigar_element_t** const cigar_buffer_sentinel,
    const cigar_t cigar_element_type,const uint64_t element_length);
void matches_cigar_buffer_add_mismatch(
    cigar_element_t** const cigar_buffer_sentinel,const uint8_t mismatch);

/*
 * CIGAR Vector Handling
 */
void matches_cigar_vector_append_insertion(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,
    const uint64_t indel_length,const cigar_attr_t attributes);
void matches_cigar_vector_append_deletion(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,
    const uint64_t indel_length,const cigar_attr_t attributes);
void matches_cigar_vector_append_match(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,
    const uint64_t match_length,const cigar_attr_t attributes);
void matches_cigar_vector_append_mismatch(
    vector_t* const cigar_vector,uint64_t* const current_cigar_length,
    const uint8_t mismatch,const cigar_attr_t attributes);
void matches_cigar_vector_append_cigar_element(
    vector_t* const cigar_vector,uint64_t* const cigar_length,cigar_element_t* const cigar_element);

/*
 * CIGAR Vector Utils
 */
void matches_cigar_reverse(
    vector_t* const cigar_vector,const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
void matches_cigar_reverse_colorspace(
    vector_t* const cigar_vector,const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);

uint64_t matches_cigar_compute_event_distance(
    const vector_t* const cigar_vector,const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
uint64_t matches_cigar_compute_edit_distance(
    const vector_t* const cigar_vector,const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
uint64_t matches_cigar_compute_edit_distance__excluding_clipping(
    const vector_t* const cigar_vector,const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
uint64_t matches_cigar_compute_matching_bases(
    const vector_t* const cigar_vector,const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);

int64_t matches_cigar_element_effective_length(const cigar_element_t* const cigar_element);
int64_t matches_cigar_effective_length(
    vector_t* const cigar_vector,const uint64_t cigar_offset,
    const uint64_t cigar_length);

/*
 * CIGAR Vector Compare
 */
int matches_cigar_cmp(
    vector_t* const cigar_vector_match0,match_trace_t* const match0,
    vector_t* const cigar_vector_match1,match_trace_t* const match1);

/*
 * Display
 */
void match_cigar_print(
    FILE* const stream,vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length);

/*
 * Error Messages
 */
#define GEM_ERROR_MATCHES_CIGAR_ZERO_LENGTH "Matches. CIGAR length cannot be zero"

#endif /* MATCHES_CIGAR_H_ */
