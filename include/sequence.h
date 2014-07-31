/*
 * PROJECT: GEMMapper
 * FILE: sequence.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple data structure to store genomic reads
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include "essentials.h"

/*
 * Checkers
 */
#define SEQUENCE_CHECK(sequence) \
  GEM_CHECK_NULL(sequence); \
  STRING_CHECK(sequence->tag); \
  STRING_CHECK(sequence->read)

#define SEQUENCE_QUALITY_IS_VALID(character) (33 <= (character))

/*
 * Sequence
 */
typedef enum { SINGLE_END, PAIRED_END1, PAIRED_END2 } sequence_end_t;
typedef struct {
  sequence_end_t end_info;
  string_t* casava_tag;
  string_t* extra_tag;
} sequence_attributes_t;
typedef struct {
  /* Sequence */
  string_t* tag;
  string_t* read;
  string_t* qualities;
  /* Attributes */
  sequence_attributes_t attributes;
} sequence_t;

/*
 * Constructor
 */
GEM_INLINE sequence_t* sequence_new(void);
GEM_INLINE void sequence_clear(sequence_t* const sequence);
GEM_INLINE void sequence_delete(sequence_t* const sequence);

/*
 * Accessors
 */
GEM_INLINE uint64_t sequence_get_length(sequence_t* const sequence);

GEM_INLINE void sequence_set_tag(sequence_t* const sequence,char* const text,const uint64_t length);
GEM_INLINE char* sequence_get_tag(sequence_t* const sequence);

GEM_INLINE void sequence_set_read(sequence_t* const sequence,char* const text,const uint64_t length);
GEM_INLINE char* sequence_get_read(sequence_t* const sequence);

GEM_INLINE void sequence_set_qualities(sequence_t* const sequence,char* const text,const uint64_t length);
GEM_INLINE char* sequence_get_qualities(sequence_t* const sequence);

GEM_INLINE bool sequence_has_qualities(const sequence_t* const sequence);
GEM_INLINE bool sequence_has_casava_tag(const sequence_t* const sequence);
GEM_INLINE bool sequence_has_extra_tag(const sequence_t* const sequence);

GEM_INLINE sequence_end_t sequence_get_end_info(const sequence_t* const sequence);

/*
 * Utils
 */
GEM_INLINE void sequence_generate_reverse(sequence_t* const sequence,sequence_t* const rev_sequence);
GEM_INLINE void sequence_generate_reverse_complement(sequence_t* const sequence,sequence_t* const rc_sequence);

#endif /* SEQUENCE_H_ */
