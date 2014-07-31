/*
 * PROJECT: GEMMapper
 * FILE: dna_string.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef DNA_STRING_H_
#define DNA_STRING_H_

#include "essentials.h"

// Range of DNA Nucleotides
#define DNA_RANGE 4
#define DNA_EXT_RANGE 7

#define DNA_RANGE_BITS 2
#define DNA_EXT_RANGE_BITS 3

/*
 * DNA Nucleotides
 */
#define DNA_CHAR_A 'A'
#define DNA_CHAR_C 'C'
#define DNA_CHAR_G 'G'
#define DNA_CHAR_T 'T'

#define DNA_CHAR_N    'N'
#define DNA_CHAR_SEP  '|'
#define DNA_CHAR_JUMP 'J'

/*
 * Encoded DNA Nucleotides
 */
#define ENC_DNA_CHAR_A 0
#define ENC_DNA_CHAR_C 1
#define ENC_DNA_CHAR_G 2
#define ENC_DNA_CHAR_T 3

#define ENC_DNA_CHAR_N    4
#define ENC_DNA_CHAR_SEP  5
#define ENC_DNA_CHAR_JUMP 6

// Orientation (strand)
typedef enum { Forward, Reverse } strand_t;

/*
 * Handy check functions
 */
extern const bool dna_table[256];
extern const bool dna_encoded_table[DNA_EXT_RANGE];
extern const bool extended_dna_table[256];
extern const bool extended_dna_encoded_table[DNA_EXT_RANGE];
extern const bool unmasked_dna_table[256];
extern const bool iupac_code_table[256];

extern const char dna_normalized_table[256];
extern const char dna_strictly_normalized_table[256];
extern const char dna_complement_table[256];
extern const uint8_t dna_encoded_complement_table[DNA_EXT_RANGE];

extern const uint8_t dna_encode_table[256];
extern const char dna_decode_table[DNA_EXT_RANGE];

extern const uint8_t dna_encoded_colorspace_table[DNA_EXT_RANGE][DNA_EXT_RANGE];

#define is_dna(character)          (dna_table[(int)(character)])
#define is_dna_encoded(character)  (dna_encoded_table[(int)(character)])
#define is_extended_dna(character)         (extended_dna_table[(int)(character)])
#define is_extended_dna_encoded(character) (extended_dna_encoded_table[(int)(character)])
#define is_unmasked_dna(character) (unmasked_dna_table[(int)(character)])
#define is_iupac_code(character)   (iupac_code_table[(int)(character)])

#define dna_normalized(character)          (dna_normalized_table[(int)(character)])
#define dna_strictly_normalized(character) (dna_strictly_normalized_table[(int)(character)])
#define dna_complement(character)          (dna_complement_table[(int)(character)])
#define dna_encoded_complement(character)  (dna_encoded_complement_table[(int)(character)])

#define dna_encode(character) (dna_encode_table[(int)(character)])
#define dna_decode(enc_char)  (dna_decode_table[(int)(enc_char)])

#define dna_encoded_colorspace(enc_char_0,enc_char_1) (dna_encoded_colorspace_table[(int)(enc_char_0)][(int)(enc_char_1)])

///*
// * Reverse-Complement functions
// */
//GEM_INLINE void dna_string_reverse_complement(dna_string* const dna_string);
//GEM_INLINE void dna_string_reverse_complement_copy(dna_string* const dna_string_dst,dna_string* const dna_string_src);
//GEM_INLINE dna_string* dna_string_reverse_complement_dup(dna_string* const dna_string);

#endif /* DNA_STRING_H_ */
