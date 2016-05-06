/*
 * PROJECT: GEMMapper
 * FILE: sequence.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple data structure to store genomic reads
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include "utils/essentials.h"

/*
 * Checkers
 */
#define SEQUENCE_QUALITY_IS_VALID(character) (33 <= (character))

/*
 * Sequence
 */
typedef enum { single_end=0, paired_end1=1, paired_end2=2 } sequence_end_t;
typedef struct {
  /* Sequence */
  string_t tag;
  string_t read;
  string_t qualities;
  /* Attributes */
  bool has_qualities;
  sequence_end_t end_info;
  string_t casava_tag;
  string_t extra_tag;
} sequence_t;

/*
 * Constructor
 */
void sequence_init(sequence_t* const sequence);
void sequence_init_mm(sequence_t* const sequence,mm_stack_t* const mm_stack);
void sequence_clear(sequence_t* const sequence);
void sequence_destroy(sequence_t* const sequence);

/*
 * Accessors
 */
uint64_t sequence_get_length(sequence_t* const sequence);

char* sequence_get_tag(sequence_t* const sequence);
void sequence_set_tag(
    sequence_t* const sequence,
    char* const text,
    const uint64_t length);

char* sequence_get_read(sequence_t* const sequence);
void sequence_set_read(
    sequence_t* const sequence,
    char* const text,
    const uint64_t length);

char* sequence_get_qualities(sequence_t* const sequence);
void sequence_set_qualities(
    sequence_t* const sequence,
    char* const text,
    const uint64_t length);

bool sequence_has_qualities(const sequence_t* const sequence);
bool sequence_has_casava_tag(const sequence_t* const sequence);
bool sequence_has_extra_tag(const sequence_t* const sequence);

sequence_end_t sequence_get_end_info(const sequence_t* const sequence);

/*
 * Utils
 */
bool sequence_equals(sequence_t* const sequence_a,sequence_t* const sequence_b);
void sequence_generate_reverse(sequence_t* const sequence,sequence_t* const rev_sequence);
void sequence_generate_reverse_complement(sequence_t* const sequence,sequence_t* const rc_sequence);

/*
 * Display
 */
void sequence_print(FILE* const stream,sequence_t* const sequence);

#endif /* SEQUENCE_H_ */
