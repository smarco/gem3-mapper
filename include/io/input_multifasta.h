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

#ifndef INPUT_MULTIFASTA_H_
#define INPUT_MULTIFASTA_H_

#include "utils/essentials.h"
#include "archive/locator.h"
#include "io/input_fasta.h"

/*
 * Utils
 */
#define MFASTA_VALID_CHARACTER(character) \
  ((character)>=BANG && (character)<=TILDE)

#define MFASTA_VALID_INITIAL_CHARACTER(character) \
  (MFASTA_VALID_CHARACTER(character) && (character)!=STAR && (character)!=EQUAL)

/*
 * MultiFASTA Parsing State
 */
typedef enum {
  Expecting_sequence,
  Expecting_tag,
  Reading_sequence
} multifasta_read_state_t;

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
 * Multi-FASTA File
 */
typedef struct {
  /* Multi-FASTA File */
  char* file_name;                        // File name
  fm_t* file_manager;
  uint64_t line_no;                       // Current line number
  /* Parsing State */
  input_multifasta_state_t parsing_state; // Parsing state
} input_multifasta_file_t;

/*
 * MultiFASTA Setup
 */
void input_multifasta_file_open(
    input_multifasta_file_t* const input_multifasta_file,
    char* const input_multifasta_file_name);
void input_multifasta_file_close(
    input_multifasta_file_t* const input_multifasta_file);

/*
 * MultiFASTA parsing state
 */
void input_multifasta_state_clear(
    input_multifasta_state_t* const parsing_state);
void input_multifasta_state_reset_interval(
    input_multifasta_state_t* const parsing_state);
void input_multifasta_state_begin_sequence(
    input_multifasta_state_t* const parsing_state);
uint64_t input_multifasta_get_text_sequence_length(
    input_multifasta_state_t* const parsing_state);

/*
 * Errors
 */
#define GEM_ERROR_MULTIFASTA_EMPTY "MultiFASTA '%s' is empty"
#define GEM_ERROR_MULTIFASTA_READING "MultiFASTA parsing (%s:%"PRIu64"). Error reading stream"
#define GEM_ERROR_MULTIFASTA_BEGINNING_TAG "MultiFASTA parsing (%s:%"PRIu64"). must start with a tag '>'"
#define GEM_ERROR_MULTIFASTA_TAG_EMPTY "MultiFASTA parsing (%s:%"PRIu64"). Tag is empty"
#define GEM_ERROR_MULTIFASTA_SEQ_EMPTY "MultiFASTA parsing (%s:%"PRIu64"). Expecting sequence (Empty sequence)"
#define GEM_ERROR_MULTIFASTA_INVALID_CHAR "MultiFASTA parsing (%s:%"PRIu64"). Invalid character found ('%c')"

#endif /* INPUT_MULTIFASTA_H_ */
