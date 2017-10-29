/*
 * PROJECT: GEM-Tools library
 * FILE: gt_segmented_sequence.c
 * DATE: 3/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 *
 */

#include "gt_segmented_sequence.h"

#define GT_SEQ_ARCHIVE_NUM_BLOCKS 15000
#define GT_SEQ_ARCHIVE_BLOCK_SIZE GT_BUFFER_SIZE_256K

/*
 * SegmentedSEQ Constructor
 */
GT_INLINE gt_segmented_sequence* gt_segmented_sequence_new(void) {
  gt_segmented_sequence* sequence = gt_alloc(gt_segmented_sequence);
  sequence->blocks = gt_vector_new(GT_SEQ_ARCHIVE_NUM_BLOCKS,sizeof(gt_compact_dna_string*));
  sequence->sequence_total_length = 0;
  sequence->seq_name = gt_string_new(10);
  return sequence;
}
GT_INLINE void gt_segmented_sequence_clear(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_VECTOR_ITERATE(sequence->blocks,block,block_num,uint64_t*) {
    if (*block) gt_free(*block);
  }
  gt_vector_clear(sequence->blocks);
  gt_string_clear(sequence->seq_name);
}
GT_INLINE void gt_segmented_sequence_delete(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_VECTOR_ITERATE(sequence->blocks,block,block_num,gt_compact_dna_string*) {
    if (*block) gt_cdna_string_delete(*block);
  }
  gt_vector_delete(sequence->blocks);
  gt_string_delete(sequence->seq_name);
  gt_free(sequence);
}
/*
 * SegmentedSEQ handler
 */
GT_INLINE void gt_segmented_sequence_set_name(gt_segmented_sequence* const sequence,char* const seq_name,const uint64_t seq_name_length) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_string_set_nstring(sequence->seq_name,seq_name,seq_name_length);
}
GT_INLINE char* gt_segmented_sequence_get_name(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  return gt_string_get_string(sequence->seq_name);
}

GT_INLINE gt_compact_dna_string* gt_segmented_sequence_get_block(gt_segmented_sequence* const sequence,const uint64_t position) {
  const uint64_t num_block = position/GT_SEQ_ARCHIVE_BLOCK_SIZE;
  const uint64_t blocks_used = gt_vector_get_used(sequence->blocks);
  // Allocate new blocks (if needed)
  if (num_block>=blocks_used) {
    uint64_t i;
    for (i=blocks_used;i<num_block;++i) {
      gt_vector_insert(sequence->blocks,NULL,gt_compact_dna_string*);
    }
    gt_compact_dna_string* const block = gt_cdna_string_new(GT_SEQ_ARCHIVE_BLOCK_SIZE);
    gt_vector_insert(sequence->blocks,block,gt_compact_dna_string*);
    return block;
  } else {
    gt_compact_dna_string* block = *gt_vector_get_elm(sequence->blocks,num_block,gt_compact_dna_string*);
    if (!block) {
      block = gt_cdna_string_new(GT_SEQ_ARCHIVE_BLOCK_SIZE);
      gt_vector_set_elm(sequence->blocks,num_block,gt_compact_dna_string*,block);
    }
    return block;
  }
}
GT_INLINE char gt_segmented_sequence_get_char_at(gt_segmented_sequence* const sequence,const uint64_t position) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_SEGMENTED_SEQ_POSITION_CHECK(sequence,position);
  const uint64_t pos_in_block = position%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  return gt_cdna_string_get_char_at(gt_segmented_sequence_get_block(sequence,position),pos_in_block);
}
GT_INLINE void gt_segmented_sequence_set_char_at(gt_segmented_sequence* const sequence,const uint64_t position,const char character) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  const uint64_t pos_in_block = position%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  // Adjust sequence total length
  if (position>=sequence->sequence_total_length) sequence->sequence_total_length = position+1;
  // Set character in compact dna string
  gt_cdna_string_set_char_at(gt_segmented_sequence_get_block(sequence,position),pos_in_block,character);
}
GT_INLINE void gt_segmented_sequence_append_string(gt_segmented_sequence* const sequence,const char* const string,const uint64_t length) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  uint64_t block_free_space = GT_SEQ_ARCHIVE_BLOCK_SIZE-(sequence->sequence_total_length%GT_SEQ_ARCHIVE_BLOCK_SIZE);
  uint64_t current_length = sequence->sequence_total_length;
  uint64_t chars_written = 0;
  while (chars_written < length) {
    const uint64_t chunk_size = ((length-chars_written)<block_free_space) ?
        length-chars_written : block_free_space;
    gt_cdna_string_append_string(gt_segmented_sequence_get_block(sequence,current_length),string+chars_written,chunk_size);
    chars_written += chunk_size;
    current_length += chunk_size;
    block_free_space = GT_SEQ_ARCHIVE_BLOCK_SIZE;
  }
  sequence->sequence_total_length = current_length;
}

GT_INLINE gt_status gt_segmented_sequence_get_sequence(
    gt_segmented_sequence* const sequence,const uint64_t position,const uint64_t length,gt_string* const string) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_SEGMENTED_SEQ_POSITION_CHECK(sequence,position);
  GT_STRING_CHECK(string);
  GT_ZERO_CHECK(length);
  // Clear string
  gt_string_clear(string);
  // Check position
  if (gt_expect_false(position >= sequence->sequence_total_length)) return GT_SEQUENCE_POS_OUT_OF_RANGE;
  // Retrieve String
  uint64_t i=0;
  gt_segmented_sequence_iterator sequence_iterator;
  gt_segmented_sequence_new_iterator(sequence,position,GT_ST_FORWARD,&sequence_iterator);
  while (i<length && !gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
    gt_string_append_char(string,gt_segmented_sequence_iterator_next(&sequence_iterator));
    ++i;
  }
  gt_string_append_eos(string);
  return (i==length) ? GT_SEQUENCE_OK : GT_SEQUENCE_CHUNK_OUT_OF_RANGE;
}

/*
 * SegmentedSEQ Iterator
 */
GT_INLINE void gt_segmented_sequence_new_iterator(
    gt_segmented_sequence* const sequence,const uint64_t position,gt_string_traversal const direction,
    gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_NULL_CHECK(sequence_iterator);
  // Set iterator
  sequence_iterator->sequence = sequence;
  if (gt_expect_true(position<sequence->sequence_total_length)) {
    gt_segmented_sequence_iterator_seek(sequence_iterator,position,direction);
  } else {
    sequence_iterator->global_pos = position; // EOS
  }
}
GT_INLINE void gt_segmented_sequence_iterator_seek(
    gt_segmented_sequence_iterator* const sequence_iterator,const uint64_t position,gt_string_traversal const direction) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  GT_SEGMENTED_SEQ_POSITION_CHECK(sequence_iterator->sequence,position);
  // Set iterator direction
  sequence_iterator->direction = direction;
  // Set sequence location fields
  sequence_iterator->global_pos = position;
  sequence_iterator->local_pos = position%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  // Init sequence locator
  gt_compact_dna_string* const cdna_string = gt_segmented_sequence_get_block(sequence_iterator->sequence,position);
  gt_cdna_string_new_iterator(cdna_string,sequence_iterator->local_pos,direction,&sequence_iterator->cdna_string_iterator);
}
GT_INLINE bool gt_segmented_sequence_iterator_eos(gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  return (sequence_iterator->global_pos>=sequence_iterator->sequence->sequence_total_length);
}
GT_INLINE char gt_segmented_sequence_iterator_following(gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  // Seek to proper position (load block if necessary)
  if (sequence_iterator->local_pos==GT_SEQ_ARCHIVE_BLOCK_SIZE) {
    gt_segmented_sequence_iterator_seek(sequence_iterator,sequence_iterator->global_pos,sequence_iterator->direction);
    gt_check(sequence_iterator->global_pos%GT_SEQ_ARCHIVE_BLOCK_SIZE != 0,ALG_INCONSISNTENCY);
  }
  // Update position
  ++sequence_iterator->global_pos;
  ++sequence_iterator->local_pos;
  // Return next!
  return gt_cdna_string_iterator_next(&sequence_iterator->cdna_string_iterator);
}
GT_INLINE char gt_segmented_sequence_iterator_previous(gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  // Seek to proper position (load block if necessary)
  if (sequence_iterator->local_pos==-1) {
    gt_segmented_sequence_iterator_seek(sequence_iterator,sequence_iterator->global_pos,sequence_iterator->direction);
    gt_check(sequence_iterator->local_pos%GT_SEQ_ARCHIVE_BLOCK_SIZE ==
        sequence_iterator->global_pos%GT_SEQ_ARCHIVE_BLOCK_SIZE,ALG_INCONSISNTENCY);
  }
  // Update position
  --sequence_iterator->global_pos;
  --sequence_iterator->local_pos;
  // Return next!
  return gt_cdna_string_iterator_next(&sequence_iterator->cdna_string_iterator);
}
GT_INLINE char gt_segmented_sequence_iterator_next(gt_segmented_sequence_iterator* const sequence_iterator) {
  return (sequence_iterator->direction==GT_ST_FORWARD) ?
    gt_segmented_sequence_iterator_following(sequence_iterator) :
    gt_segmented_sequence_iterator_previous(sequence_iterator);
}
