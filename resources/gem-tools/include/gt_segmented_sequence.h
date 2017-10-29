/*
 * PROJECT: GEM-Tools library
 * FILE: gt_segmented_sequence.h
 * DATE: 3/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   // TODO
 */

#ifndef GT_SEGMENTED_SEQUENCE_H_
#define GT_SEGMENTED_SEQUENCE_H_

#include "gt_essentials.h"

#include "gt_dna_string.h"
#include "gt_compact_dna_string.h"

/*
 * Codes gt_status
 */
#define GT_SEQUENCE_OK 0

#define GT_SEQUENCE_POS_OUT_OF_RANGE 11
#define GT_SEQUENCE_CHUNK_OUT_OF_RANGE 12

/*
 * Data types
 */
typedef struct {
  gt_string* seq_name;
  gt_vector* blocks; /* (gt_compact_dna_string*) */
  uint64_t sequence_total_length;
} gt_segmented_sequence;

typedef struct {
  gt_segmented_sequence* sequence;
  gt_string_traversal direction;
  uint64_t global_pos;
  int64_t local_pos;
  gt_compact_dna_string_iterator cdna_string_iterator;
} gt_segmented_sequence_iterator;

/*
 * Checkers
 */
#define GT_SEGMENTED_SEQ_CHECK(segmented_sequence) \
  GT_NULL_CHECK(segmented_sequence); \
  GT_NULL_CHECK(segmented_sequence->blocks); \
  GT_STRING_CHECK(segmented_sequence->seq_name)
#define GT_SEGMENTED_SEQ_POSITION_CHECK(segmented_sequence,position) \
  gt_fatal_check(position>=segmented_sequence->sequence_total_length, \
      SEGMENTED_SEQ_IDX_OUT_OF_RANGE,position,segmented_sequence->sequence_total_length);
#define GT_SEGMENTED_SEQ_ITERATOR_CHECK(segmented_sequence_iterator) \
  GT_SEGMENTED_SEQ_CHECK(segmented_sequence_iterator->sequence)

/*
 * SegmentedSEQ Constructor
 */
GT_INLINE gt_segmented_sequence* gt_segmented_sequence_new(void);
GT_INLINE void gt_segmented_sequence_clear(gt_segmented_sequence* const sequence);
GT_INLINE void gt_segmented_sequence_delete(gt_segmented_sequence* const sequence);
/*
 * SegmentedSEQ Sequence handler
 */
GT_INLINE void gt_segmented_sequence_set_name(gt_segmented_sequence* const sequence,char* const seq_name,const uint64_t seq_name_length);
GT_INLINE char* gt_segmented_sequence_get_name(gt_segmented_sequence* const sequence);

GT_INLINE char gt_segmented_sequence_get_char_at(gt_segmented_sequence* const sequence,const uint64_t position);
GT_INLINE void gt_segmented_sequence_set_char_at(gt_segmented_sequence* const sequence,const uint64_t position,const char character);
GT_INLINE void gt_segmented_sequence_append_string(gt_segmented_sequence* const sequence,const char* const string,const uint64_t length);

GT_INLINE gt_status gt_segmented_sequence_get_sequence(
    gt_segmented_sequence* const sequence,const uint64_t position,const uint64_t length,gt_string* const string);
/*
 * SegmentedSEQ Iterator
 */
GT_INLINE void gt_segmented_sequence_new_iterator(
    gt_segmented_sequence* const sequence,const uint64_t position,gt_string_traversal const direction,
    gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE void gt_segmented_sequence_iterator_seek(
    gt_segmented_sequence_iterator* const sequence_iterator,const uint64_t position,gt_string_traversal const direction);
GT_INLINE bool gt_segmented_sequence_iterator_eos(gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE char gt_segmented_sequence_iterator_next(gt_segmented_sequence_iterator* const sequence_iterator);

/*
 * Error Messages
 */
#define GT_ERROR_SEGMENTED_SEQ_IDX_OUT_OF_RANGE "Error accessing segmented sequence. Index %"PRIu64" out out range [0,%"PRIu64")"

#endif /* GT_SEGMENTED_SEQUENCE_H_ */
