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
typedef enum { cigar_match=0, cigar_mismatch=1, cigar_ins=2, cigar_del=3, cigar_soft_trim=4, cigar_null=5 } cigar_t;
typedef struct {
  int32_t indel_length;   // Indel length
  uint8_t* indel_text;    // Indel-text (Reference text in case of insertion, NULL otherwise)
} cigar_indel_t;
typedef struct {
  cigar_t type;           // Match, Mismatch, insertion or deletion
  union {
    int32_t match_length; // Match length
    uint8_t mismatch;     // Mismatch base
    cigar_indel_t indel;  // Pointer to the insert
  };
} cigar_element_t;
/*
 * Interval Match
 */
typedef struct {
  /* Index */
  uint64_t lo;            // SA Lo-Position
  uint64_t hi;            // SA Hi-Position
  /* Sequence */
  uint8_t* text;          // Pointer to the matching-text
  uint64_t length;        // Length of the matching-text
  uint64_t distance;      // Distance of the alignment
  int32_t swg_score;      // SWG Distance (score) of the alignment
  strand_t strand;        // Mapping Strand (make a difference between matches found by searching F/R)
} match_interval_t;
/*
 * Trace Match
 */
typedef struct {
  // Error degree of the region matching
  uint64_t error; // TODO
  // Coordinates of the region
  uint64_t key_begin;
  uint64_t key_end;
  uint64_t text_begin;
  uint64_t text_end;
} region_matching_t;
typedef struct {
  /* Match */
  char* sequence_name;     // Sequence name (After decoding)
  uint64_t index_position; // Position of the match in the index. Global index => (seq1|seq2|...|)
  uint64_t text_position;  // Position of the match in the text. Local text => seq15
  uint64_t trace_offset;   // Chained list of matching blocks
  uint64_t distance;       // Distance of the alignment
  int32_t swg_score;       // SWG Distance (score) of the alignment
  strand_t strand;         // Mapping Strand
  /* Score */
  uint8_t mapq_score;      // MAPQ Score
  /* CIGAR */
  uint64_t cigar_buffer_offset;
  uint64_t cigar_length;
  int64_t effective_length;
#ifdef GEM_DEBUG
  uint64_t num_regions_matching;
  region_matching_t* regions_matching;
#endif
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
  /* Search-matches state */
  uint64_t max_complete_stratum;
  /* Text Collection Buffer */
  text_collection_t* text_collection;       // Stores text-traces (candidates/matches/regions/...)
  /* Matches Counters */
  vector_t* counters;                       // Global counters
  uint64_t min_counter_value;
  uint64_t max_counter_value;
  /* Interval Matches */
  vector_t* interval_matches;               // (match_interval_t)
  /* Global Position Matches */
  vector_t* global_matches;                 // Global Matches (match_trace_t)
  ihash_t* begin_gmatches;                  // Begin position (of the aligned match) in the text-space
  ihash_t* end_gmatches;                    // End position (of the aligned match) in the text-space
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
} matches_t;

/*
 * Setup
 */
GEM_INLINE matches_t* matches_new();
GEM_INLINE void matches_configure(matches_t* const matches,text_collection_t* const text_collection);
GEM_INLINE void matches_clear(matches_t* const matches);
GEM_INLINE void matches_delete(matches_t* const matches);

/*
 * Counters
 */
GEM_INLINE uint64_t matches_counters_get_min(matches_t* const matches);
GEM_INLINE uint64_t matches_counters_get_max(matches_t* const matches);
GEM_INLINE uint64_t matches_counters_compact(matches_t* const matches);

/*
 * Trace-Matches
 */
GEM_INLINE cigar_element_t* match_trace_get_cigar_array(const matches_t* const matches,const match_trace_t* const match_trace);
GEM_INLINE uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace);
GEM_INLINE uint64_t match_trace_get_distance(const match_trace_t* const match_trace);
GEM_INLINE int64_t match_trace_get_effective_length(
    matches_t* const matches,const uint64_t read_length,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length);

/*
 * Adding Matches
 */
GEM_INLINE uint64_t* matches_lookup_match(
    matches_t* const matches,const uint64_t begin_position,const uint64_t effective_length);
GEM_INLINE void matches_add_match_trace_t(
    matches_t* const matches,match_trace_t* const match_trace,
    const bool update_counters,mm_stack_t* const mm_stack);

GEM_INLINE void matches_add_interval_match(
    matches_t* const matches,const uint64_t lo,const uint64_t hi,
    const uint64_t length,const uint64_t distance,const strand_t strand);
GEM_INLINE void matches_add_interval_set(
    matches_t* const matches,interval_set_t* const interval_set,
    const uint64_t length,const strand_t strand);

GEM_INLINE void matches_hint_allocate_match_trace(matches_t* const matches,const uint64_t num_matches_trace_to_add);
GEM_INLINE void matches_hint_allocate_match_interval(matches_t* const matches,const uint64_t num_matches_interval_to_add);

/*
 * CIGAR Handling
 */
GEM_INLINE void matches_cigar_buffer_add_cigar_element(
    cigar_element_t** const cigar_buffer_sentinel,const cigar_t cigar_element_type,
    const uint64_t element_length,uint8_t* const indel_text);

GEM_INLINE void matches_cigar_buffer_append_indel(
    vector_t* const cigar_buffer,uint64_t* const current_cigar_length,
    const cigar_t cigar_element_type,const uint64_t element_length,uint8_t* const indel_text);
GEM_INLINE void matches_cigar_buffer_append_match(
    vector_t* const cigar_buffer,uint64_t* const current_cigar_length,
    const uint64_t match_length);
GEM_INLINE void matches_cigar_buffer_append_mismatch(
    vector_t* const cigar_buffer,uint64_t* const current_cigar_length,
    const cigar_t cigar_element_type,const uint8_t mismatch);

GEM_INLINE void matches_cigar_reverse(
    matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length);
GEM_INLINE void matches_cigar_reverse_colorspace(
    matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length);

GEM_INLINE uint64_t matches_cigar_calculate_edit_distance(
    const matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length);
GEM_INLINE uint64_t matches_cigar_calculate_edit_distance__excluding_clipping(
    const matches_t* const matches,const uint64_t cigar_buffer_offset,const uint64_t cigar_length);

/*
 * Status
 */
GEM_INLINE bool matches_is_mapped(const matches_t* const matches);

/*
 * Sorting Matches
 */
GEM_INLINE void matches_sort_by_distance(matches_t* const matches);
GEM_INLINE void matches_sort_by_mapq_score(matches_t* const matches);
GEM_INLINE void matches_sort_by_sequence_name__position(matches_t* const matches);

/*
 * Restore Point
 */
GEM_INLINE void matches_restore_point_save(matches_t* const matches);
GEM_INLINE void matches_restore_point_rollback(matches_t* const matches);

/*
 * Error Messages
 */
#define GEM_ERROR_MATCHES_CIGAR_ZERO_LENGTH "Matches. CIGAR length cannot be zero"

//
//#define MATCHES_GET_KEY_INFO(POS_MATCH,MATCHES) ((keys_info*)vector_get_mem(MATCHES->qbuf_keys_info)+POS_MATCH->key_id)
//#define FMI_MATCHES_GET_NUM_MATCHES(matches) (vector_get_used(matches->rbuf_pos) + matches->num_int_matches)
//#define MATCHES_IS_DUPLICATED(mA_pos,mB_pos,mA_len,mB_len,mA_key_id,mB_key_id) ((mA_key_id==mB_key_id) && ((mA_pos)==(mB_pos) || (mA_pos)+(mA_len)==(mB_pos)+(mB_len)))


#endif /* MATCHES_H_ */
