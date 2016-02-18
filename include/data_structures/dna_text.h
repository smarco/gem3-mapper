/*
 * PROJECT: GEMMapper
 * FILE: compact_dna_text.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Provides functionality to handle a compact representation of a 8-chars alphabet text
 */

#ifndef DNA_TEXT_H_
#define DNA_TEXT_H_

#include "utils/essentials.h"

/*
 * Range of DNA Nucleotides
 */
#define DNA_RANGE     4
#define DNA__N_RANGE  5
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

/*
 * Translation tables
 */
extern const bool dna_canonical_table[256];
extern const bool dna_canonical_encoded_table[DNA_EXT_RANGE];
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

extern const char dna_bisulfite_C2T_table[256];
extern const char dna_bisulfite_G2A_table[256];
extern const uint8_t dna_encoded_bisulfite_C2T_table[DNA_EXT_RANGE];
extern const uint8_t dna_encoded_bisulfite_G2A_table[DNA_EXT_RANGE];

extern const uint8_t dna_encoded_colorspace_table[DNA_EXT_RANGE][DNA_EXT_RANGE];

/*
 * Translation functions
 */
#define is_dna_canonical(character)          (dna_canonical_table[(int)(character)])
#define is_dna_canonical_encoded(character)  (dna_canonical_encoded_table[(int)(character)])
#define is_dna(character)                    (dna_table[(int)(character)])
#define is_dna_encoded(character)            (dna_encoded_table[(int)(character)])
#define is_extended_dna(character)           (extended_dna_table[(int)(character)])
#define is_extended_dna_encoded(enc_char)    (extended_dna_encoded_table[(int)(enc_char)])
#define is_unmasked_dna(character)           (unmasked_dna_table[(int)(character)])
#define is_iupac_code(character)             (iupac_code_table[(int)(character)])

#define dna_normalized(character)            (dna_normalized_table[(int)(character)])
#define dna_strictly_normalized(character)   (dna_strictly_normalized_table[(int)(character)])
#define dna_complement(character)            (dna_complement_table[(int)(character)])
#define dna_encoded_complement(enc_char)     (dna_encoded_complement_table[(int)(enc_char)])

#define dna_bisulfite_C2T(character)         (dna_bisulfite_C2T_table[(int)(character)])
#define dna_bisulfite_G2A(character)         (dna_bisulfite_G2A_table[(int)(character)])
#define dna_encoded_bisulfite_C2T(enc_char)  (dna_encoded_bisulfite_C2T_table[(int)(enc_char)])
#define dna_encoded_bisulfite_G2A(enc_char)  (dna_encoded_bisulfite_G2A_table[(int)(enc_char)])

#define dna_encode(character)                (dna_encode_table[(int)(character)])
#define dna_decode(enc_char)                 (dna_decode_table[(int)(enc_char)])

#define dna_encoded_colorspace(enc_char_0,enc_char_1) (dna_encoded_colorspace_table[(int)(enc_char_0)][(int)(enc_char_1)])

/*
 * Orientation (strand)
 */
typedef enum { Forward=0, Reverse=UINT64_MAX } strand_t;             // DNA strand type
typedef enum {
  bs_strand_none = 0,
  bs_strand_C2T = 1,
  bs_strand_G2A = 2,
  bs_strand_mixed = UINT64_MAX
} bs_strand_t; // Bisulfite strand type

/*
 * DNA-Text
 */
typedef struct _dna_text_t dna_text_t;

/*
 * Setup/Loader
 */
dna_text_t* dna_text_read_mem(mm_t* const memory_manager);
void dna_text_delete(dna_text_t* const dna_text);

/*
 * Builder
 */
dna_text_t* dna_text_new(const uint64_t dna_text_length);
dna_text_t* dna_text_padded_new(
    const uint64_t dna_text_length,
    const uint64_t init_padding,
    const uint64_t end_padding);
void dna_text_write_chunk(
    fm_t* const output_file_manager,
    dna_text_t* const dna_text,
    const uint64_t chunk_length);
void dna_text_write(
    fm_t* const output_file_manager,
    dna_text_t* const dna_text);

/*
 * Accessors
 */
uint64_t dna_text_get_length(const dna_text_t* const dna_text);
void dna_text_set_length(
    dna_text_t* const dna_text,
    const uint64_t length);
uint64_t dna_text_get_size(const dna_text_t* const dna_text);
uint8_t dna_text_get_char(
    const dna_text_t* const dna_text,
    const uint64_t position);
void dna_text_set_char(
    const dna_text_t* const dna_text,
    const uint64_t position,
    const uint8_t enc_char);
uint8_t* dna_text_get_text(const dna_text_t* const dna_text);
uint8_t* dna_text_retrieve_sequence(
    const dna_text_t* const dna_text,
    const uint64_t position,
    const uint64_t length,
    mm_stack_t* const mm_stack);

/*
 * Utils
 */
strand_t dna_text_strand_get_complement(const strand_t strand);

/*
 * Display
 */
void dna_text_print(
    FILE* const stream,
    dna_text_t* const dna_text,
    const uint64_t length);
void dna_text_print_content(
    FILE* const stream,
    dna_text_t* const dna_text);
void dna_text_pretty_print_content(
    FILE* const stream,
    dna_text_t* const dna_text,
    const uint64_t width);

// Buffer Display
void dna_buffer_print(
    FILE* const stream,
    const uint8_t* const dna_buffer,
    const uint64_t dna_buffer_length,
    const bool print_reverse);

/*
 * Errors
 */
#define GEM_ERROR_DNA_TEXT_WRONG_MODEL_NO "DNA-text. Wrong DNA-text Model %"PRIu64" (Expected model %"PRIu64")"
#define GEM_ERROR_DNA_TEXT_OOR "DNA-text. Requested position (%"PRIu64") out of range [0,%"PRIu64")"

#endif /* DNA_TEXT_H_ */
