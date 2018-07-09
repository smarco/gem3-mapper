/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Input module allows parsing of MultiFASTA files
 */

#include "io/input_multifasta.h"
#include "text/dna_text.h"

/*
 * MultiFASTA Setup
 */
void input_multifasta_file_open(
    input_multifasta_file_t* const input_multifasta_file,
    char* const input_multifasta_file_name) {
  input_multifasta_file->file_name = input_multifasta_file_name;
  input_multifasta_file->file_manager = fm_open_file(input_multifasta_file_name,FM_READ);
  input_multifasta_file->line_no = 0;
  input_multifasta_state_clear(&input_multifasta_file->parsing_state);
}
void input_multifasta_file_close(
    input_multifasta_file_t* const input_multifasta_file) {
  fm_close(input_multifasta_file->file_manager);
}
/*
 * MultiFASTA parsing state
 */
void input_multifasta_state_clear(
    input_multifasta_state_t* const parsing_state) {
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
void input_multifasta_state_reset_interval(
    input_multifasta_state_t* const parsing_state) {
  /* Text */
  parsing_state->text_sequence_length += parsing_state->text_interval_length;
  parsing_state->text_interval_length = 0;
  /* Index */
  parsing_state->index_sequence_length += parsing_state->index_interval_length;
  parsing_state->index_interval_length = 0;
}
void input_multifasta_state_begin_sequence(
    input_multifasta_state_t* const parsing_state) {
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
uint64_t input_multifasta_get_text_sequence_length(
    input_multifasta_state_t* const parsing_state) {
  return (parsing_state->text_sequence_length+parsing_state->text_interval_length+parsing_state->ns_pending);
}

