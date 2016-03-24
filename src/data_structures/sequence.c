/*
 * PROJECT: GEMMapper
 * FILE: sequence.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple data structure to store genomic reads
 */

#include "data_structures/sequence.h"
#include "data_structures/dna_text.h"

#define SEQUENCE_TAG_INITIAL_LENGTH 80
#define SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH 40
#define SEQUENCE_INITIAL_LENGTH 200

/*
 * Constructor
 */
void sequence_init(sequence_t* const sequence) {
  string_init(&sequence->tag,SEQUENCE_TAG_INITIAL_LENGTH);
  string_init(&sequence->read,SEQUENCE_INITIAL_LENGTH);
  string_init(&sequence->qualities,SEQUENCE_INITIAL_LENGTH);
  sequence->end_info = single_end;
  string_init(&sequence->casava_tag,SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH);
  string_init(&sequence->extra_tag,SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH);
}
void sequence_init_mm(sequence_t* const sequence,mm_stack_t* const mm_stack) {
  string_init_mm(&sequence->tag,SEQUENCE_TAG_INITIAL_LENGTH,mm_stack);
  string_init_mm(&sequence->read,SEQUENCE_INITIAL_LENGTH,mm_stack);
  string_init_mm(&sequence->qualities,SEQUENCE_INITIAL_LENGTH,mm_stack);
  sequence->end_info = single_end;
  string_init_mm(&sequence->casava_tag,SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH,mm_stack);
  string_init_mm(&sequence->extra_tag,SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH,mm_stack);
}
void sequence_clear(sequence_t* const sequence) {
  string_clear(&sequence->tag);
  string_clear(&sequence->read);
  string_clear(&sequence->qualities);
  sequence->end_info = single_end;
  string_clear(&sequence->casava_tag);
  string_clear(&sequence->extra_tag);
}
void sequence_destroy(sequence_t* const sequence) {
  string_destroy(&sequence->tag);
  string_destroy(&sequence->read);
  string_destroy(&sequence->qualities);
  string_destroy(&sequence->casava_tag);
  string_destroy(&sequence->extra_tag);
}
/*
 * Accessors
 */
uint64_t sequence_get_length(sequence_t* const sequence) {
  return string_get_length(&sequence->read);
}
char* sequence_get_tag(sequence_t* const sequence) {
  return string_get_buffer(&sequence->tag);
}
void sequence_set_tag(
    sequence_t* sequence,
    char* const text,
    const uint64_t length) {
  string_set_buffer(&sequence->tag,text,length);
}
char* sequence_get_read(sequence_t* const sequence) {
  return string_get_buffer(&sequence->read);
}
void sequence_set_read(
    sequence_t* sequence,
    char* const text,
    const uint64_t length) {
  string_set_buffer(&sequence->read,text,length);
}
char* sequence_get_qualities(sequence_t* const sequence) {
  return string_get_buffer(&sequence->qualities);
}
void sequence_set_qualities(
    sequence_t* sequence,
    char* const text,
    const uint64_t length) {
  string_set_buffer(&sequence->qualities,text,length);
}
bool sequence_has_qualities(const sequence_t* const sequence) {
  return !string_is_null(&sequence->qualities);
}
bool sequence_has_casava_tag(const sequence_t* const sequence) {
  return !string_is_null(&sequence->casava_tag);
}
bool sequence_has_extra_tag(const sequence_t* const sequence) {
  return !string_is_null(&sequence->extra_tag);
}
sequence_end_t sequence_get_end_info(const sequence_t* const sequence) {
  return sequence->end_info;
}
/*
 * Utils
 */
bool sequence_equals(sequence_t* const sequence_a,sequence_t* const sequence_b) {
  return string_equals(&sequence_a->read,&sequence_b->read);
}
void sequence_generate_reverse(sequence_t* const sequence,sequence_t* const rev_sequence) {
  // Prepare rc_string (Read)
  const uint64_t seq_buffer_length = string_get_length(&sequence->read);
  string_resize(&rev_sequence->read,seq_buffer_length);
  string_clear(&rev_sequence->read);
  // Reverse Read
  int64_t pos;
  const char* const seq_buffer = string_get_buffer(&sequence->read);
  for (pos=seq_buffer_length-1;pos>=0;--pos) {
    string_append_char(&rev_sequence->read,seq_buffer[pos]);
  }
  string_append_eos(&rev_sequence->read);
  // Reverse Qualities
  if (sequence_has_qualities(sequence)) {
    string_copy_reverse(&rev_sequence->qualities,&sequence->qualities);
  }
}
void sequence_generate_reverse_complement(sequence_t* const sequence,sequence_t* const rc_sequence) {
  // Prepare rc_string (Read)
  const uint64_t seq_buffer_length = string_get_length(&sequence->read);
  string_resize(&rc_sequence->read,seq_buffer_length);
  string_clear(&rc_sequence->read);
  // Reverse-Complement Read
  int64_t pos;
  const char* const seq_buffer = string_get_buffer(&sequence->read);
  for (pos=seq_buffer_length-1;pos>=0;--pos) {
    string_append_char(&rc_sequence->read,dna_complement(seq_buffer[pos]));
  }
  string_append_eos(&rc_sequence->read);
  // Reverse Qualities
  if (sequence_has_qualities(sequence)) {
    string_copy_reverse(&rc_sequence->qualities,&sequence->qualities);
  }
}
/*
 * Display
 */
void sequence_print(FILE* const stream,sequence_t* const sequence) {
  tab_fprintf(stream,"[GEM]>Sequence\n");
  if (sequence_has_qualities(sequence)) {
    fprintf(stream,"@%s\n",sequence_get_tag(sequence));
    fprintf(stream,"%s\n+\n",sequence_get_read(sequence));
    fprintf(stream,"%s\n",sequence_get_qualities(sequence));
  } else {
    fprintf(stream,">%s\n",sequence_get_tag(sequence));
    fprintf(stream,"%s\n",sequence_get_read(sequence));
  }
}

