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

#include "commons.h"
#include "fm.h"
#include "mm.h"
#include "mm_stack.h"

/*
 * DNA-Text Model
 */
#define DNA_TEXT_RAW
// #define DNA_TEXT_COMPACT

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

/*
 * Checkers
 */
#define DNA_TEXT_CHECK(dna_text) GEM_CHECK_NULL(dna_text)

/*
 * Orientation (strand)
 */
typedef enum { Forward, Reverse } strand_t;
/*
 * DNA-Text
 */
typedef struct _dna_text_t dna_text_t;
typedef struct _dna_text_builder_t dna_text_builder_t;
///*
// * DNA-Text Stats
// */
//typedef struct {
//  /* Nucleotides */
//  stats_vector_t* nucleotides;
//  stats_vector_t* dimers;
//  stats_vector_t* trimers;
//  stats_vector_t* tetramers;
//  /* Runs */
//  stats_vector_t* runs;
//  /* Distance between nucleotide repetitions */
//  uint64_t* last_position;
//  stats_vector_t* distance_nucleotides;
//  /* Abundance => How many different nucleotides are there in each @SIZE-bucket (@SIZE chars) */
//  stats_vector_t* abundance_256nt;
//  stats_vector_t* abundance_128nt;
//  stats_vector_t* abundance_64nt;
//  stats_vector_t* abundance_32nt;
//  stats_vector_t* abundance_16nt;
//  stats_vector_t* abundance_8nt;
//  /* Internals */ // TODO
//} dna_text_stats_t;

/*
 * Setup/Loader
 */
GEM_INLINE dna_text_t* dna_text_read(fm_t* const file_manager);
GEM_INLINE dna_text_t* dna_text_read_mem(mm_t* const memory_manager);
GEM_INLINE void dna_text_delete(dna_text_t* const dna_text);

/*
 * Builder
 */
GEM_INLINE dna_text_builder_t* dna_text_builder_new(const uint64_t dna_text_length);
GEM_INLINE dna_text_builder_t* dna_text_builder_padded_new(
    const uint64_t dna_text_length,const uint64_t init_padding,const uint64_t end_padding);
GEM_INLINE void dna_text_builder_write(
    fm_t* const output_file_manager,dna_text_builder_t* const dna_text);
GEM_INLINE void dna_text_builder_delete(dna_text_builder_t* const dna_text);

GEM_INLINE uint64_t dna_text_builder_get_length(const dna_text_builder_t* const dna_text);
GEM_INLINE void dna_text_builder_set_length(dna_text_builder_t* const dna_text,const uint64_t length);
GEM_INLINE uint8_t dna_text_builder_get_char(const dna_text_builder_t* const dna_text,const uint64_t position);
GEM_INLINE void dna_text_builder_set_char(const dna_text_builder_t* const dna_text,const uint64_t position,const uint8_t enc_char);
GEM_INLINE uint8_t* dna_text_builder_get_buffer(const dna_text_builder_t* const dna_text);

/*
 * Accessors
 */
GEM_INLINE uint64_t dna_text_get_length(const dna_text_t* const dna_text);
GEM_INLINE uint8_t* dna_text_get_buffer(const dna_text_t* const dna_text);
GEM_INLINE uint8_t* dna_text_retrieve_sequence(
    const dna_text_t* const dna_text,const uint64_t position,const uint64_t length,
    mm_stack_t* const mm_stack);

/*
 * Display
 */
GEM_INLINE void dna_text_print(FILE* const stream,dna_text_t* const dna_text);
GEM_INLINE void dna_text_builder_print(FILE* const stream,dna_text_builder_t* const dna_text);
GEM_INLINE void dna_text_builder_print_content(FILE* const stream,dna_text_builder_t* const dna_text);
GEM_INLINE void dna_text_builder_pretty_print_content(FILE* const stream,dna_text_builder_t* const dna_text,const uint64_t width);

///*
// * DNA Text [Stats]
// */
//GEM_INLINE dna_text_stats_t* dna_text_stats_new();
//GEM_INLINE void dna_text_stats_delete(dna_text_stats_t* const dna_text_stats);
//
//// Calculate Stats
//GEM_INLINE void dna_text_stats_record(dna_text_stats_t* const dna_text_stats,const uint8_t char_enc);
//
//// Display
//GEM_INLINE void dna_text_stats_print(FILE* const stream,dna_text_stats_t* const dna_text_stats);

/*
 * Errors
 */
#define GEM_ERROR_DNA_TEXT_OOR "DNA-text. Requested position (%lu) out of range [0,%lu)"

#endif /* DNA_TEXT_H_ */
