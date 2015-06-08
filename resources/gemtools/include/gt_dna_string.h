/*
 * PROJECT: GEM-Tools library
 * FILE: gt_dna_string.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_DNA_STRING_H_
#define GT_DNA_STRING_H_

#include "gt_essentials.h"

#define GT_DNA_RANGE 5

#define GT_DNA_CHAR_A 'A'
#define GT_DNA_CHAR_C 'C'
#define GT_DNA_CHAR_G 'G'
#define GT_DNA_CHAR_T 'T'
#define GT_DNA_CHAR_N 'N'

#define gt_dna_string gt_string

// Orientation (strand)
typedef enum { FORWARD, REVERSE, UNKNOWN } gt_strand;

typedef struct {
  gt_dna_string* dna_string;
  uint64_t current_pos;
} gt_dna_string_iterator;

typedef enum { GT_ST_FORWARD, GT_ST_BACKWARD } gt_string_traversal;

/*
 * Handy check functions
 */
extern const bool gt_dna[256];
extern const bool gt_iupac_code[256];

extern const char gt_dna_normalized[256];
extern const char gt_dna_strictly_normalized[256];
extern const char gt_complement_table[256];

#define gt_is_dna(character)         (gt_dna[(int)(character)])
#define gt_is_iupac_code(character)  (gt_iupac_code[(int)(character)])

#define gt_get_dna_normalized(character) (gt_dna_normalized[(int)(character)])
#define gt_get_dna_strictly_normalized(character) (gt_dna_strictly_normalized[(int)(character)])
#define gt_get_complement(character) (gt_complement_table[(int)(character)])

/*
 * Checkers
 */
#define GT_DNA_STRING_CHECK(dna_string) \
  GT_NULL_CHECK(dna_string); \
  GT_STRING_CHECK(dna_string)
#define GT_DNA_STRING_ITERATOR_CHECK(dna_string_iterator) \
  GT_NULL_CHECK(dna_string_iterator); \
  GT_DNA_STRING_CHECK(dna_string_iterator->dna_string)

/*
 * Constructor
 */
#define gt_dna_string_new gt_string_new
#define gt_dna_string_resize gt_string_resize
#define gt_dna_string_clear gt_string_clear
#define gt_dna_string_delete gt_string_delete

/*
 * DNA String handler
 */
GT_INLINE bool gt_dna_string_is_dna_string(gt_dna_string* const dna_string);

GT_INLINE char gt_dna_string_get_char_at(gt_dna_string* const dna_string,const uint64_t pos); // TODO
GT_INLINE void gt_dna_string_set_char_at(gt_dna_string* const dna_string,const uint64_t pos,const char character); // TODO

#define gt_dna_string_get_string gt_string_get_string
#define gt_dna_string_get_length gt_string_get_length

GT_INLINE void gt_dna_string_set_string(gt_dna_string* const dna_string,char* const dna_string_src);
GT_INLINE void gt_dna_string_set_nstring(gt_dna_string* const dna_string,char* const dna_string_src,const uint64_t length);

GT_INLINE void gt_dna_string_reverse_complement(gt_dna_string* const dna_string);
GT_INLINE void gt_dna_string_reverse_complement_copy(gt_dna_string* const dna_string_dst,gt_dna_string* const dna_string_src);
GT_INLINE gt_dna_string* gt_dna_string_reverse_complement_dup(gt_dna_string* const dna_string);

/*
 * DNA String Iterator
 */
GT_INLINE void gt_dna_string_new_iterator(
    gt_dna_string* const dna_string,const uint64_t pos,gt_string_traversal const direction,
    gt_dna_string_iterator* const dna_string_iterator);
GT_INLINE void gt_dna_string_iterator_seek(gt_dna_string_iterator* const dna_string_iterator,const uint64_t pos);
GT_INLINE bool gt_dna_string_iterator_eos(gt_dna_string_iterator* const dna_string_iterator);
GT_INLINE char gt_dna_string_iterator_next(gt_dna_string_iterator* const dna_string_iterator);
GT_INLINE char gt_dna_string_iterator_previous(gt_dna_string_iterator* const dna_string_iterator);

#endif /* GT_DNA_STRING_H_ */
