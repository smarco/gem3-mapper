/*
 * PROJECT: GEMMapper
 * FILE: input_multifasta_parser.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for MULTIFASTA files
 */

#include "io/input_multifasta_parser.h"
#include "data_structures/dna_text.h"

/*
 * MultiFASTA parsing state
 */
void input_multifasta_state_clear(input_multifasta_state_t* const parsing_state) {
  /* Parsing State */
  parsing_state->multifasta_read_state = Expecting_tag;
  /* Sequence components */
  parsing_state->tag_id = 0;                 // Current sequence TAG-ID
  parsing_state->text_position = 0;          // Current position of the current sequence (from MultiFASTA)
  parsing_state->index_position = 0;         // Current position of generated index
  /* Text */
  parsing_state->ns_pending = 0;             // Accumulated Ns
  parsing_state->text_interval_length = 0;   // Length of the text-interval
  parsing_state->text_sequence_length = 0;   // Length of the text-sequence (Sum of intervals)
  /* Index */
  parsing_state->index_interval_length = 0;  // Length of the index-interval
  parsing_state->index_sequence_length = 0;  // Length of the index-sequence (Sum of intervals)
  parsing_state->last_char = DNA_CHAR_SEP;   // Last character printed in the index
  parsing_state->interval_type = locator_interval_chromosomal_assembly; // Current interval type
  parsing_state->strand = Forward;           // Current strand
  parsing_state->bs_strand = bs_strand_none; // Current BS-Strand
}
void input_multifasta_state_reset_interval(input_multifasta_state_t* const parsing_state) {
  /* Text */
  parsing_state->text_sequence_length += parsing_state->text_interval_length;
  parsing_state->text_interval_length = 0;
  /* Index */
  parsing_state->index_sequence_length += parsing_state->index_interval_length;
  parsing_state->index_interval_length = 0;
}
void input_multifasta_state_begin_sequence(input_multifasta_state_t* const parsing_state) {
  parsing_state->multifasta_read_state = Expecting_sequence;
  /* Text */
  parsing_state->text_position = 0;
  parsing_state->ns_pending = 0;
  parsing_state->text_interval_length = 0;
  parsing_state->text_sequence_length = 0;
  /* Index */
  parsing_state->index_interval_length = 0;
  parsing_state->index_sequence_length = 0;
}
uint64_t input_multifasta_get_text_sequence_length(input_multifasta_state_t* const parsing_state) {
  return (parsing_state->text_sequence_length+parsing_state->text_interval_length+parsing_state->ns_pending);
}
/*
 * MultiFASTA file parsing
 */
void input_multifasta_parse_tag(input_file_t* const input_multifasta,string_t* const tag) {
  // Prepare to add a new TAG
  string_clear(tag);
  // Check empty tag
  gem_cond_fatal_error(!input_file_next_char(input_multifasta),
      MULTIFASTA_TAG_EMPTY,PRI_input_file_content(input_multifasta));
  // Read tag
  char current_char = input_file_get_current_char(input_multifasta);
  while (!IS_ANY_EOL(current_char) &&
         !MFASTA_IS_ANY_TAG_SEPARATOR(current_char) &&
         input_file_next_char(input_multifasta)) {
    string_append_char(tag,current_char); // Append to tag
    current_char = input_file_get_current_char(input_multifasta); // Next Character
  }
  if (!IS_ANY_EOL(current_char) && !MFASTA_IS_ANY_TAG_SEPARATOR(current_char)) {
    string_append_char(tag,current_char);
  }
  // Check empty and append EOS
  gem_cond_fatal_error(string_get_length(tag)==0,
      MULTIFASTA_TAG_EMPTY,PRI_input_file_content(input_multifasta)); // Empty TAG
  string_append_eos(tag);
  // Skip the rest of the line
  while (!IS_ANY_EOL(current_char) && input_file_next_char(input_multifasta)) {
    current_char = input_file_get_current_char(input_multifasta);
  }
  input_file_skip_eol(input_multifasta);
}
void input_multifasta_skip_tag(input_file_t* const input_multifasta) {
  // Check empty tag
  gem_cond_fatal_error(!input_file_next_char(input_multifasta),
      MULTIFASTA_TAG_EMPTY,PRI_input_file_content(input_multifasta));
  // Read tag
  char current_char = input_file_get_current_char(input_multifasta);
  while (!IS_ANY_EOL(current_char) &&
         !MFASTA_IS_ANY_TAG_SEPARATOR(current_char) &&
         input_file_next_char(input_multifasta)) {
    current_char = input_file_get_current_char(input_multifasta); // Next Character
  }
  // Skip the rest of the line
  while (!IS_ANY_EOL(current_char) && input_file_next_char(input_multifasta)) {
    current_char = input_file_get_current_char(input_multifasta);
  }
  input_file_skip_eol(input_multifasta);
}

