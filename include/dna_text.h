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

#include "essentials.h"
#include "dna_string.h"

#include "stats_matrix.h"
#include "stats_vector.h"

/*
 * Checkers
 */
#define DNA_TEXT_CHECK(dna_text) \
  GEM_CHECK_NULL(dna_text)

/*
 *#define DNA_TEXT_CHECK(dna_text) \
 *  GEM_CHECK_NULL(dna_text); \
 *  GEM_CHECK_NULL(dna_text->text)
 */

/*
 * DNA Text
 */
typedef struct _dna_text_t dna_text_t;
/*
 * DNA Text Stats
 */
typedef struct {
  /* Nucleotides */
  stats_vector_t* nucleotides;
  stats_vector_t* dimers;
  stats_vector_t* trimers;
  stats_vector_t* tetramers;
  /* Runs */
  stats_vector_t* runs;
  /* Distance between nucleotide repetitions */
  uint64_t* last_position;
  stats_vector_t* distance_nucleotides;
  /* Abundance => How many different nucleotides are there in each @SIZE-bucket (@SIZE chars) */
  stats_vector_t* abundance_256nt;
  stats_vector_t* abundance_128nt;
  stats_vector_t* abundance_64nt;
  stats_vector_t* abundance_32nt;
  stats_vector_t* abundance_16nt;
  stats_vector_t* abundance_8nt;
  /* Internals */
  // TODO
} dna_text_stats_t;

/*
 * DNA Text
 */
GEM_INLINE dna_text_t* dna_text_new(const uint64_t dna_text_length);
GEM_INLINE dna_text_t* dna_text_padded_new(
    const uint64_t dna_text_length,const uint64_t init_padding,const uint64_t end_padding);
GEM_INLINE void dna_text_delete(dna_text_t* const dna_text);

GEM_INLINE dna_text_t* dna_text_read(fm_t* const file_manager);
GEM_INLINE dna_text_t* dna_text_read_mem(mm_t* const memory_manager);

GEM_INLINE void dna_text_write(fm_t* const output_file_manager,dna_text_t* const dna_text);

// Accessors
GEM_INLINE uint64_t dna_text_get_length(const dna_text_t* const dna_text);
GEM_INLINE void dna_text_set_length(dna_text_t* const dna_text,const uint64_t length);
GEM_INLINE uint8_t* dna_text_get_buffer(const dna_text_t* const dna_text);
GEM_INLINE uint8_t dna_text_get_char(const dna_text_t* const dna_text,const uint64_t position);
GEM_INLINE void dna_text_set_char(const dna_text_t* const dna_text,const uint64_t position,const uint8_t enc_char);

// Display
GEM_INLINE void dna_text_print(FILE* const stream,dna_text_t* const dna_text);
GEM_INLINE void dna_text_print_content(FILE* const stream,dna_text_t* const dna_text);
GEM_INLINE void dna_text_pretty_print_content(FILE* const stream,dna_text_t* const dna_text,const uint64_t width);

/*
 * DNA Text [Stats]
 */
GEM_INLINE dna_text_stats_t* dna_text_stats_new();
GEM_INLINE void dna_text_stats_delete(dna_text_stats_t* const dna_text_stats);

// Calculate Stats
GEM_INLINE void dna_text_stats_record(dna_text_stats_t* const dna_text_stats,const uint8_t char_enc);

// Display
GEM_INLINE void dna_text_stats_print(FILE* const stream,dna_text_stats_t* const dna_text_stats);

/*
 * Errors
 */
#define GEM_ERROR_DNA_TEXT_OOR "dna-text. Requested position (%lu) out of range [0,%lu)"

#endif /* DNA_TEXT_H_ */
