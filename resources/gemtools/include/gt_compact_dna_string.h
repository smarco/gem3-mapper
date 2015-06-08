/*
 * PROJECT: GEM-Tools library
 * FILE: gt_compact_dna_string.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Bitmap (compact) representation of DNA-strings (using 8 characters alphabet)
 */

#ifndef GT_COMPACT_DNA_STRING_H_
#define GT_COMPACT_DNA_STRING_H_

#include "gt_commons.h"
#include "gt_dna_string.h"

/*
 * Blocks encoding dimensions
 */
#define GT_CDNA_BLOCK_CHARS 64
#define GT_CDNA_BLOCK_BITMAP_SIZE 8
#define GT_CDNA_BLOCK_BITMAPS 3
#define GT_CDNA_BLOCK_SIZE (GT_CDNA_BLOCK_BITMAP_SIZE*GT_CDNA_BLOCK_BITMAPS)

/*
 * CDNA Bitmaps Encoding
 */
#define GT_CDNA_ENC_CHAR_A 0
#define GT_CDNA_ENC_CHAR_C 1
#define GT_CDNA_ENC_CHAR_G 2
#define GT_CDNA_ENC_CHAR_T 3
#define GT_CDNA_ENC_CHAR_N 4

#define GT_CDNA_ZERO_MASK       0xFFFFFFFFFFFFFFFEull
#define GT_CDNA_ONE_MASK        0x0000000000000001ull
#define GT_CDNA_ONE_LAST_MASK   0x8000000000000000ull
#define GT_CDNA_ZERO_LAST_MASK  0x7FFFFFFFFFFFFFFFull

#define GT_CDNA_EXTRACT_MASK      GT_CDNA_ONE_MASK
#define GT_CDNA_EXTRACT_LAST_MASK GT_CDNA_ONE_LAST_MASK

/*
 * CDNA DataStructures
 */
typedef struct {
  uint64_t* bitmaps;
  uint64_t allocated;
  uint64_t length;
} gt_compact_dna_string;

typedef struct {
  gt_compact_dna_string* cdna_string;
  uint64_t current_pos;
  uint64_t current_pos_mod;
  gt_string_traversal direction;
  uint64_t* current_bitmap;
  uint64_t bm_0;
  uint64_t bm_1;
  uint64_t bm_2;
} gt_compact_dna_string_iterator;

extern const char gt_cdna_decode[8];
extern const uint8_t gt_cdna_encode[256];

#define gt_cdna_decode(enc_char)  gt_cdna_decode[enc_char]
#define gt_cdna_encode(character) gt_cdna_encode[(uint8_t)character]

/*
 * Checkers
 */
#define GT_COMPACT_DNA_STRING_CHECK(cdna_string) \
  GT_NULL_CHECK(cdna_string); \
  GT_NULL_CHECK(cdna_string->bitmaps)
#define GT_COMPACT_DNA_STRING_POSITION_CHECK(cdna_string,position) \
  gt_check(position>=cdna_string->length,CDNA_IT_OUT_OF_RANGE,position,cdna_string->length);
#define GT_COMPACT_DNA_STRING_ITERATOR_CHECK(cdna_string_iterator) \
  GT_NULL_CHECK(cdna_string_iterator); \
  GT_COMPACT_DNA_STRING_CHECK(cdna_string_iterator->cdna_string); \
  GT_NULL_CHECK(cdna_string_iterator->current_bitmap)

/*
 * Constructor
 */
GT_INLINE gt_compact_dna_string* gt_cdna_string_new(const uint64_t initial_chars);
GT_INLINE void gt_cdna_string_resize(gt_compact_dna_string* const cdna_string,const uint64_t num_chars);
GT_INLINE void gt_cdna_string_clear(gt_compact_dna_string* const cdna_string);
GT_INLINE void gt_cdna_string_delete(gt_compact_dna_string* const cdna_string);

/*
 * Handlers
 */
GT_INLINE char gt_cdna_string_get_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos);
GT_INLINE void gt_cdna_string_set_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos,const char character);
GT_INLINE uint64_t gt_cdna_string_get_length(gt_compact_dna_string* const cdna_string);
GT_INLINE void gt_cdna_string_append_string(gt_compact_dna_string* const cdna_string,const char* const string,const uint64_t length);

/*
 * Compact DNA String Sequence Iterator
 */
GT_INLINE void gt_cdna_string_new_iterator(
    gt_compact_dna_string* const cdna_string,const uint64_t position,gt_string_traversal const direction,
    gt_compact_dna_string_iterator* const cdna_string_iterator);
GT_INLINE void gt_cdna_string_iterator_seek(
    gt_compact_dna_string_iterator* const cdna_string_iterator,
    const uint64_t pos,gt_string_traversal const direction);
GT_INLINE bool gt_cdna_string_iterator_eos(gt_compact_dna_string_iterator* const cdna_string_iterator);
GT_INLINE char gt_cdna_string_iterator_next(gt_compact_dna_string_iterator* const cdna_string_iterator);

/*
 * Error Messages
 */
#define GT_ERROR_CDNA_IT_OUT_OF_RANGE "Error seeking sequence. Index %"PRIu64" out out range [0,%"PRIu64")"

#endif /* GT_COMPACT_DNA_STRING_H_ */
