/*
 * PROJECT: GEMMapper
 * FILE: matches.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure to store alignment matches {sequence,position,strand,...}
 */

#ifndef MATCHES_H_
#define MATCHES_H_

#include "essentials.h"
#include "dna_string.h"
#include "interval_set.h"
#include "text_collection.h"

/*
 * Checkers
 */
#define MATCHES_CHECK(matches) GEM_CHECK_NULL(matches)
#define MATCH_TRACE_CHECK(match) GEM_CHECK_NULL(match)

/*
 * CIGAR (Mismatches/Indels)
 */
typedef enum { cigar_match, cigar_mismatch, cigar_ins, cigar_del, cigar_soft_trim } cigar_t;
typedef struct {
  cigar_t type;    // Match, Mismatch, insertion or deletion
  int32_t length; // Match length, Mismatch base or indel length
} cigar_element_t;

/*
 * Interval Match
 */
typedef struct {
  /* Index */
  uint64_t lo;           // SA Lo-Position
  uint64_t hi;           // SA Hi-Position
  /* Sequence */
  char* text;            // Pointer to the matching-text
  uint64_t length;       // Length of the matching-text
  uint64_t distance;     // Levenshtein score of the alignment
  strand_t strand;       // Mapping Strand (make a difference between matches found by searching F/R)
} match_interval_t;
/*
 * Trace Match
 */
typedef struct {
  /* Match */
  char* sequence_name;   // Sequence name (After decoding)
  uint64_t position;     // Position of the match
  uint64_t trace_offset; // Chained list of matching blocks
  uint64_t distance;     // Edit distance of the alignment
  strand_t strand;       // Mapping Strand
  /* Score */
  int64_t score;
  /* CIGAR */
  uint64_t cigar_buffer_offset;
  uint64_t cigar_length;
} match_trace_t;

///*
// * Local Match // TODO
// */
//typedef struct _match_trace_local_t match_trace_local_t;
//struct _match_trace_local_t {
//  /* Match */
//  uint64_t trace_offset;           // Chained list of matching blocks
//  uint64_t distance;               // Edit distance of the alignment
//  strand_t strand;                 // Mapping Strand
//  match_trace_local_t* next_block; // Next block     (Chained list of matching local blocks)
//  match_trace_local_t* prev_block; // Previous block (Chained list of matching local blocks)
//  uint64_t sequence_begin;         // Begin Position of the sequence block matching
//  uint64_t sequence_end;           // End Position of the sequence block matching
//  /* CIGAR */
//  uint64_t cigar_buffer_offset;
//  uint64_t cigar_length;
//};


/*
 * Matches
 */
typedef struct {
  /* Text Collection Buffer */
  text_collection_t* text_collection;       // Stores text-traces (candidates/matches/regions/...)
  /* Matches Counters */
  vector_t* counters;                       // Global counters
  /* Interval Matches */
  vector_t* interval_matches;               // (match_interval_t)
  /* Global Position Matches */
  vector_t* global_matches;                 // Global Matches (match_trace_t)
  ihash_t* start_gmatches;                  // Starting position in the Text-space
//  /* Local Position Matches */ // TODO
//  vector_t* local_matches;               // (match_trace_t)
//  ihash_t* start_position_gmatches;      // Starting position in the Text-space
//  ihash_t* end_position_gmatches;        // End position in the Text-space
//  vector_t* local_matches_block_buffer;  // Match Position Block Buffer (match_position_block_t)
  /* Buffers */
  vector_t* cigar_buffer;                  // CIGAR operations buffer (cigar_element_t)
  /* Restore Point (RP) */
  vector_t* rp_counters;
  uint64_t rp_interval_matches_used;
  uint64_t rp_global_matches_used;
  // uint64_t rp_local_matches_used; // TODO
  uint64_t rp_cigar_buffer_used;
  /* Pre-computed Data */
  uint64_t total_matches;
  uint64_t last_computed_interval_matches_used;
  uint64_t last_computed_global_matches_used;
} matches_t;

/*
 * Setup
 */
GEM_INLINE matches_t* matches_new();
GEM_INLINE void matches_clear(matches_t* const matches);
GEM_INLINE void matches_delete(matches_t* const matches);

/*
 * Counters
 */
GEM_INLINE uint64_t matches_get_num_matches(matches_t* const matches);
GEM_INLINE uint64_t matches_counters_compact(matches_t* const matches);
GEM_INLINE uint64_t matches_counters_get_min_matching_stratum(matches_t* const matches);

/*
 * Trace-Matches
 */
GEM_INLINE cigar_element_t* match_trace_get_cigar_array(const matches_t* const matches,const match_trace_t* const match_trace);
GEM_INLINE uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace);

/*
 * Adding Matches
 */
GEM_INLINE match_trace_t* matches_lookup_match(matches_t* const matches,const uint64_t position);
GEM_INLINE void matches_add_match_trace_t(
    matches_t* const matches,match_trace_t* const match_trace);
GEM_INLINE void matches_add_match_trace_(
    matches_t* const matches,const uint64_t trace_offset,
    const uint64_t position,const uint64_t distance,const strand_t strand);

GEM_INLINE void matches_add_interval_match(
    matches_t* const matches,
    const uint64_t hi,const uint64_t lo,
    const uint64_t length,const uint64_t distance,const strand_t strand);
GEM_INLINE void matches_add_interval_set(
    matches_t* const matches,interval_set_t* const interval_set);

/*
 * Handling matches
 */
GEM_INLINE void matches_reverse_CIGAR(
    matches_t* const matches,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length);
GEM_INLINE void matches_reverse_CIGAR_colorspace(
    matches_t* const matches,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length);
GEM_INLINE uint64_t matches_get_effective_length(
    matches_t* const matches,const uint64_t read_length,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length);

/*
 * Sorting Matches
 */
GEM_INLINE void matches_sort_by_distance(matches_t* const matches);

/*
 * Restore Point
 */
GEM_INLINE void matches_restore_point_save(matches_t* const matches);
GEM_INLINE void matches_restore_point_rollback(matches_t* const matches);

/*
 * Error Messages
 */
// #define GEM_ERROR_MATCHES_ ""

//
//#define MATCHES_GET_KEY_INFO(POS_MATCH,MATCHES) ((keys_info*)vector_get_mem(MATCHES->qbuf_keys_info)+POS_MATCH->key_id)
//#define FMI_MATCHES_GET_NUM_MATCHES(matches) (vector_get_used(matches->rbuf_pos) + matches->num_int_matches)
//#define MATCHES_IS_DUPLICATED(mA_pos,mB_pos,mA_len,mB_len,mA_key_id,mB_key_id) ((mA_key_id==mB_key_id) && ((mA_pos)==(mB_pos) || (mA_pos)+(mA_len)==(mB_pos)+(mB_len)))


#endif /* MATCHES_H_ */
