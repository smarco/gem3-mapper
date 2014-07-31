/*
 * PROJECT: GEMMapper
 * FILE: sequence.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple data structure to store genomic reads
 */

#include "sequence.h"
#include "dna_string.h"

#define SEQUENCE_TAG_INITIAL_LENGTH 80
#define SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH 40
#define SEQUENCE_INITIAL_LENGTH 200

/*
 * Constructor
 */
GEM_INLINE sequence_t* sequence_new(void) {
  sequence_t* sequence = mm_alloc(sequence_t);
  sequence->tag = string_new(SEQUENCE_TAG_INITIAL_LENGTH);
  sequence->read = string_new(SEQUENCE_INITIAL_LENGTH);
  sequence->qualities = string_new(SEQUENCE_INITIAL_LENGTH);
  sequence->attributes.end_info = SINGLE_END;
  sequence->attributes.casava_tag = string_new(SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH);
  sequence->attributes.extra_tag = string_new(SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH);
  return sequence;
}
GEM_INLINE void sequence_clear(sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  string_clear(sequence->tag);
  string_clear(sequence->read);
  string_clear(sequence->qualities);
  sequence->attributes.end_info = SINGLE_END;
  string_clear(sequence->attributes.casava_tag);
  string_clear(sequence->attributes.extra_tag);
}
GEM_INLINE void sequence_delete(sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  string_delete(sequence->tag);
  string_delete(sequence->read);
  string_delete(sequence->qualities);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t sequence_get_length(sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return string_get_length(sequence->read);
}
GEM_INLINE void sequence_set_tag(sequence_t* sequence,char* const text,const uint64_t length) {
  SEQUENCE_CHECK(sequence);
  string_set_buffer(sequence->tag,text,length);
}
GEM_INLINE char* sequence_get_tag(sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return string_get_buffer(sequence->tag);
}
GEM_INLINE void sequence_set_read(sequence_t* sequence,char* const text,const uint64_t length) {
  SEQUENCE_CHECK(sequence);
  string_set_buffer(sequence->read,text,length);
}
GEM_INLINE char* sequence_get_read(sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return string_get_buffer(sequence->read);
}
GEM_INLINE void sequence_set_qualities(sequence_t* sequence,char* const text,const uint64_t length) {
  SEQUENCE_CHECK(sequence);
  string_set_buffer(sequence->qualities,text,length);
}
GEM_INLINE char* sequence_get_qualities(sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return string_get_buffer(sequence->qualities);
}
GEM_INLINE bool sequence_has_qualities(const sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return !string_is_null(sequence->qualities);
}
GEM_INLINE bool sequence_has_casava_tag(const sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return !string_is_null(sequence->attributes.casava_tag);
}
GEM_INLINE bool sequence_has_extra_tag(const sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return !string_is_null(sequence->attributes.extra_tag);
}
GEM_INLINE sequence_end_t sequence_get_end_info(const sequence_t* const sequence) {
  SEQUENCE_CHECK(sequence);
  return sequence->attributes.end_info;
}
/*
 * Utils
 */
GEM_INLINE void sequence_generate_reverse(sequence_t* const sequence,sequence_t* const rev_sequence) {
  // Prepare rc_string (Read)
  const uint64_t seq_buffer_length = string_get_length(sequence->read);
  string_resize(rev_sequence->read,seq_buffer_length);
  string_clear(rev_sequence->read);
  // Reverse Read
  int64_t pos;
  const char* const seq_buffer = string_get_buffer(sequence->read);
  for (pos=seq_buffer_length-1;pos>=0;--pos) {
    string_append_char(rev_sequence->read,seq_buffer[pos]);
  }
  string_append_eos(rev_sequence->read);
  // Reverse Qualities
  if (sequence_has_qualities(sequence)) {
    string_copy_reverse(rev_sequence->qualities,sequence->qualities);
  }
}
GEM_INLINE void sequence_generate_reverse_complement(sequence_t* const sequence,sequence_t* const rc_sequence) {
  // Prepare rc_string (Read)
  const uint64_t seq_buffer_length = string_get_length(sequence->read);
  string_resize(rc_sequence->read,seq_buffer_length);
  string_clear(rc_sequence->read);
  // Reverse-Complement Read
  int64_t pos;
  const char* const seq_buffer = string_get_buffer(sequence->read);
  for (pos=seq_buffer_length-1;pos>=0;--pos) {
    string_append_char(rc_sequence->read,dna_complement(seq_buffer[pos]));
  }
  string_append_eos(rc_sequence->read);
  // Reverse Qualities
  if (sequence_has_qualities(sequence)) {
    string_copy_reverse(rc_sequence->qualities,sequence->qualities);
  }
}

// TODO
// TODO
// TODO
//// Reverse & Complement (not in the case of colorspace), and check if it is symmetric
// register const uint64_t key_len = search_params->key_len;
// register ch_t* const base_key = search_params->key;
// register ch_t* const base_mismatch_mask = search_params->mismatch_mask;
// register bool is_symmetric = true;
// register uint64_t i;
// vector_prepare(mpool->qbuf1, ch_t, 2*(key_len+1)); // Allocate memory for the reverse
// search_params->key = vector_get_mem(mpool->qbuf1);
// search_params->key[key_len] = 0;
// search_params->mismatch_mask = search_params->key + (key_len+1);
// search_params->mismatch_mask[key_len] = 0;
// if (archive->filter_type == Iupac_colorspace_dna) {
//   for (i = 0; i < key_len; ++i) {
//     register const uint64_t inv_pos = key_len-1-i;
//     search_params->key[i] = base_key[inv_pos];
//     is_symmetric = is_symmetric && search_params->key[i]==base_key[i];
//   }
// } else {
//   for (i = 0; i < key_len; ++i) {
//     register const uint64_t inv_pos = key_len-1-i;
//     search_params->key[i] = complement[base_key[inv_pos]];
//     is_symmetric = is_symmetric && search_params->key[i]==base_key[i];
//   }
// }
//
// if (!archive_approximate_search->is_symmetric) {
//   if (base_mismatch_mask) {
//     for (i = 0; i < key_len; ++i) {
//       search_params->mismatch_mask[i] = base_mismatch_mask[key_len-1-i];
//     }
//   } else {
//     search_params->mismatch_mask = 0;
//   }














