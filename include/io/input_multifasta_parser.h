/*
 * PROJECT: GEMMapper
 * FILE: input_multifasta_parser.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for MULTIFASTA files
 */

#ifndef INPUT_MULTIFASTA_PARSER_H_
#define INPUT_MULTIFASTA_PARSER_H_

#include "utils/essentials.h"
#include "archive/locator.h"
#include "io/input_file.h"
#include "io/input_fasta_parser.h"

/*
 * Utils
 */
#define MFASTA_IS_ANY_TAG_SEPARATOR(character) \
  ((character)==TAB || (character)==COLON || (character)==COMA || (character)==EOL || (character)==EOS)

/*
 * MultiFASTA Parsing State
 */
typedef enum { Expecting_sequence, Expecting_tag, Reading_sequence } multifasta_read_state_t;
typedef struct {
  /* Parsing State */
  multifasta_read_state_t multifasta_read_state;
  /* Sequence components */
  int64_t tag_id;                      // Current sequence TAG-ID
  uint64_t text_position;              // Current position of the current sequence (from MultiFASTA)
  uint64_t index_position;             // Current position of generated index
  /* Text */
  uint64_t ns_pending;                 // Accumulated Ns
  uint64_t text_interval_length;       // Length of the text-interval
  uint64_t text_sequence_length;       // Length of the text-sequence (Sum of intervals)
  /* Index */
  uint64_t index_interval_length;      // Length of the index-interval
  uint64_t index_sequence_length;      // Length of the index-sequence (Sum of intervals)
  char last_char;                      // Last character dumped in the index
  locator_interval_type interval_type; // Current interval type
  strand_t strand;                     // Current strand
  bs_strand_t bs_strand;               // Current BS-Strand
} input_multifasta_state_t;

/*
 * MultiFASTA parsing state
 */
void input_multifasta_state_clear(input_multifasta_state_t* const parsing_state);
void input_multifasta_state_reset_interval(input_multifasta_state_t* const parsing_state);
void input_multifasta_state_begin_sequence(input_multifasta_state_t* const parsing_state);
uint64_t input_multifasta_get_text_sequence_length(input_multifasta_state_t* const parsing_state);

/*
 * MultiFASTA file parsing
 */
void input_multifasta_parse_tag(input_file_t* const input_multifasta,string_t* const tag);
void input_multifasta_skip_tag(input_file_t* const input_multifasta);

/*
 * Errors
 */
#define GEM_ERROR_MULTIFASTA_EMPTY "MultiFASTA '%s' is empty"
#define GEM_ERROR_MULTIFASTA_BEGINNING_TAG "MultiFASTA parsing (%"PRI_input_file"). must start with a tag '>'"
#define GEM_ERROR_MULTIFASTA_TAG_EMPTY "MultiFASTA parsing (%"PRI_input_file"). Tag is empty"
#define GEM_ERROR_MULTIFASTA_SEQ_EMPTY "MultiFASTA parsing (%"PRI_input_file"). Expecting sequence (Empty sequence)"
#define GEM_ERROR_MULTIFASTA_INVALID_CHAR "MultiFASTA parsing (%"PRI_input_file"). Invalid character found ('%c')"

#endif /* INPUT_MULTIFASTA_PARSER_H_ */
