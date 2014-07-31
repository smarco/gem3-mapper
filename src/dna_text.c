/*
 * PROJECT: GEMMapper
 * FILE: dna_text.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Provides functionality to handle a compact representation of a 8-alphabet text
 */

#include "dna_text.h"

struct _dna_text_t {
  /* Text */
  uint8_t* buffer;
  uint8_t* text;
  uint64_t length;
  uint64_t allocated;
  /* MM */
  mm_t* mm_text;
  bool mm_extern;
};

/*
 * DNA Text
 */
GEM_INLINE dna_text_t* dna_text_new(const uint64_t dna_text_length) {
  dna_text_t* const dna_text = mm_alloc(dna_text_t);
  dna_text->allocated = dna_text_length;
  dna_text->buffer = mm_calloc(dna_text->allocated,uint8_t,false);
  dna_text->text = dna_text->buffer;
  dna_text->length = 0;
  dna_text->mm_text = NULL;
  dna_text->mm_extern = false;
  return dna_text;
}
GEM_INLINE dna_text_t* dna_text_padded_new(
    const uint64_t dna_text_length,const uint64_t init_padding,const uint64_t end_padding) {
  dna_text_t* const dna_text = mm_alloc(dna_text_t);
  dna_text->allocated = dna_text_length+init_padding+end_padding;
  dna_text->buffer = mm_calloc(dna_text->allocated,uint8_t,false);
  dna_text->text = dna_text->buffer + init_padding;
  dna_text->length = 0;
  dna_text->mm_text = NULL;
  dna_text->mm_extern = false;
  return dna_text;
}
GEM_INLINE void dna_text_delete(dna_text_t* const dna_text) {
  DNA_TEXT_CHECK(dna_text);
  if (dna_text->mm_text!=NULL) {
    mm_bulk_free(dna_text->mm_text);
  } else if (!dna_text->mm_extern) {
    mm_free(dna_text->buffer);
  }
  mm_free(dna_text);
}
GEM_INLINE dna_text_t* dna_text_read(fm_t* const file_manager) {
  // Alloc
  dna_text_t* const dna_text = mm_alloc(dna_text_t);
  // Read header
  dna_text->length = fm_read_uint64(file_manager);
  dna_text->allocated = dna_text->length;
  // Read Text
  dna_text->mm_extern = false;
  dna_text->mm_text = fm_load_mem(file_manager,dna_text->length*UINT8_SIZE);
  dna_text->buffer = mm_get_base_mem(dna_text->mm_text);
  dna_text->text = dna_text->buffer;
  // Return
  return dna_text;
}
GEM_INLINE dna_text_t* dna_text_read_mem(mm_t* const memory_manager) {
  // Alloc
  dna_text_t* const dna_text = mm_alloc(dna_text_t);
  // Read header
  dna_text->length = mm_read_uint64(memory_manager);
  dna_text->allocated = dna_text->length;
  // Read Text
  dna_text->mm_extern = true;
  dna_text->mm_text = NULL;
  dna_text->buffer = mm_read_mem(memory_manager,dna_text->length*UINT8_SIZE);
  dna_text->text = dna_text->buffer;
  // Return
  return dna_text;
}
GEM_INLINE void dna_text_write(fm_t* const output_file_manager,dna_text_t* const dna_text) {
  FM_CHECK(output_file_manager);
  DNA_TEXT_CHECK(dna_text);
  fm_write_uint64(output_file_manager,dna_text->length);
  fm_write_mem(output_file_manager,dna_text->text,dna_text->length*UINT8_SIZE);
}
// Accessors
GEM_INLINE uint64_t dna_text_get_length(const dna_text_t* const dna_text) {
  DNA_TEXT_CHECK(dna_text);
  return dna_text->length;
}
GEM_INLINE void dna_text_set_length(dna_text_t* const dna_text,const uint64_t length) {
  DNA_TEXT_CHECK(dna_text);
  dna_text->length = length;
}
GEM_INLINE uint8_t* dna_text_get_buffer(const dna_text_t* const dna_text) {
  DNA_TEXT_CHECK(dna_text);
  return dna_text->text;
}
GEM_INLINE uint8_t dna_text_get_char(const dna_text_t* const dna_text,const uint64_t position) {
  gem_fatal_check(position >= dna_text->allocated,DNA_TEXT_OOR,position,dna_text->allocated);
  return dna_text->text[position];
}
GEM_INLINE void dna_text_set_char(const dna_text_t* const dna_text,const uint64_t position,const uint8_t enc_char) {
  gem_fatal_check(position >= dna_text->allocated,DNA_TEXT_OOR,position,dna_text->allocated);
  dna_text->text[position] = enc_char;
}
// Display
GEM_INLINE void dna_text_print(FILE* const stream,dna_text_t* const dna_text) {
  DNA_TEXT_CHECK(dna_text);
  fprintf(stream,"[GEM]>DNA-text\n");
  fprintf(stream,"  => Architecture Raw.encoded\n");
  fprintf(stream,"  => Text.Length %lu\n",dna_text->length);
  fprintf(stream,"  => Text.Size %lu MB\n",CONVERT_B_TO_MB(dna_text->length*UINT8_SIZE));
  fflush(stream); // Flush
}
GEM_INLINE void dna_text_print_content(FILE* const stream,dna_text_t* const dna_text) {
  DNA_TEXT_CHECK(dna_text);
  const uint8_t* const enc_text = dna_text->text;
  const uint64_t text_length = dna_text->length;
  fwrite(enc_text,1,text_length,stream);
}
GEM_INLINE void dna_text_print_content_folded(FILE* const stream,dna_text_t* const dna_text,const uint64_t width) {
  DNA_TEXT_CHECK(dna_text);
  // Iterate over all indexed text
  const uint8_t* const enc_text = dna_text->text;
  const uint64_t text_length = dna_text->length;
  uint64_t i, imod=0;
  for (i=0;i<text_length;++i) {
    // Print each character
    fprintf(stream,"%c",dna_decode(enc_text[i]));
    // Print column-wise
    if (++imod==width) {
      fprintf(stream,"\n"); imod = 0;
    }
  }
  if (imod!=width) fprintf(stream,"\n");
}
/*
 * DNA Text [Stats]
 */
GEM_INLINE dna_text_stats_t* dna_text_stats_new() {
  return NULL; // TODO
}
GEM_INLINE void dna_text_stats_delete(dna_text_stats_t* const dna_text_stats) {
  // TODO
}
// Calculate Stats
GEM_INLINE void dna_text_stats_record(dna_text_stats_t* const dna_text_stats,const uint8_t char_enc) {
  // TODO
}
// Display
GEM_INLINE void dna_text_stats_print(FILE* const stream,dna_text_stats_t* const dna_text_stats) {
  // TODO
}

