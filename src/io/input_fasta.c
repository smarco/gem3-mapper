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
 *   Input module allows fast parsing of FASTA/FASTQ buffered-input files
 */

#include "io/input_fasta.h"
#include "text/dna_text.h"
#include "io/input_text.h"

/*
 * Error codes
 */
#define FASTA_ERROR_TAG_BEGINNING 20
#define FASTA_ERROR_READ_BAD_CHARACTER 30
#define FASTA_ERROR_SEPARATOR_BAD_CHARACTER 40
#define FASTA_ERROR_QUALITIES_BAD_CHARACTER 50
#define FASTA_ERROR_LENGTHS 60

/*
 * Error Messages
 */
#define GEM_ERROR_PARSE_FASTQ "Parsing FASTA/FASTQ error(%s:%"PRIu64")"
#define GEM_ERROR_PARSE_FASTQ_READ_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Input read sequence contains bad character %s"
#define GEM_ERROR_PARSE_FASTQ_QUALITIES_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Input read qualities contains invalid quality value"
#define GEM_ERROR_PARSE_FASTQ_TAG_BEGINNING "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Beginning Symbol ('>' or '@') not found. Bad syntax"
#define GEM_ERROR_PARSE_FASTQ_SEPARATOR_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Wrong separator symbol ('+'). Bad syntax"
#define GEM_ERROR_PARSE_FASTQ_LENGTHS "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Different read and qualities lengths"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * FASTQ/FASTA Tag Parser
 */
sequence_end_t input_fasta_parse_sequence_tag_chomp_end_info(string_t* const tag) {
  // Parse the end information {/1,/2}
  const char* const tag_buffer = string_get_buffer(tag);
  const uint64_t tag_length = string_get_length(tag);
  if (tag_length > 2 && tag_buffer[tag_length-2] == SLASH) {
    switch (tag_buffer[tag_length-1]) {
      case '1':
        string_set_length(tag,tag_length-2);
        return paired_end1;
      case '2':
      case '3':
        string_set_length(tag,tag_length-2);
        return paired_end2;
      default:
        break;
    }
  }
  return single_end;
}
void input_fasta_parse_sequence_tag_extra(
    sequence_t* const sequence,
    char* const tag_extra,
    const uint64_t tag_extra_length) {
  // Parse extra info
  uint64_t i = 0;
  while (i < tag_extra_length) {
    // Skip separator
    if (tag_extra[i]!=SPACE && tag_extra[i]!=TAB) break;
    ++i;
    /*
     * CASAVA 1.8 Attributes
     *   @EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACG
     */
    if (input_text_parse_count_colons_in_field(tag_extra+i)==3) {
      char* const casava_tag_begin = tag_extra+i;
      switch (casava_tag_begin[0]) {
        case '1':
          sequence->end_info = paired_end1;
          break;
        case '2':
        case '3':
          sequence->end_info = paired_end2;
          break;
        default:
          sequence->end_info = single_end;
          break;
      }
      while (i < tag_extra_length && tag_extra[i]!=SPACE && tag_extra[i]!=TAB) ++i;
      char* const casava_tag_end = tag_extra+i;
      string_set_buffer(&sequence->casava_tag,casava_tag_begin,(casava_tag_end-casava_tag_begin));
      continue; // Next!
    }
    /*
     * Extra Tag (Consider everything else as extra information appended to the tag)
     */
    char* const extra_tag_begin = tag_extra+i;
    while (i < tag_extra_length && tag_extra[i]!=SPACE && tag_extra[i]!=TAB) ++i;
    char* const extra_tag_end = tag_extra+i;
    string_set_buffer(&sequence->extra_tag,extra_tag_begin,extra_tag_end-extra_tag_begin);
    /*
     * Additional check to see if any extra attributes end in /1 /2 /3
     *     If no pair has been found yet, /1/2/3 is cut away
     *     and the tag is reset to the original tag but spaces are replaced with _
     *     E.g
     *       @SRR384920.1 HWI-ST382_0049:1:1:1217:1879/1
     *       @SRR384920.1_HWI-ST382_0049:1:1:1217:1879
     */
    const sequence_end_t end_info = input_fasta_parse_sequence_tag_chomp_end_info(&sequence->extra_tag);
    if (end_info==paired_end1 || end_info==paired_end2) {
      // Append to the tag an '_' plus the extra information
      string_append_char(&sequence->tag,UNDERSCORE);
      string_append_string(&sequence->tag,&sequence->extra_tag);
      // Replace all spaces
      const uint64_t tag_length = string_get_length(&sequence->tag);
      char* const tag_buffer = string_get_buffer(&sequence->tag);
      for (i=0;i<tag_length;i++) {
        if (tag_buffer[i]==SPACE) tag_buffer[i] = UNDERSCORE;
      }
      string_append_eos(&sequence->tag);
      string_clear(&sequence->extra_tag);
      sequence->end_info = end_info; // Set pair info
    }
  }
}
int input_fasta_parse_sequence_tag(
    buffered_input_file_t* const buffered_input,
    sequence_t* const sequence) {
  // Read TAG line
  string_t* const tag_line = &sequence->extra_tag;
  buffered_input_file_get_line(buffered_input,tag_line);
  // Parse beginning character
  const uint64_t line_length = string_get_length(tag_line)-1; // No EOL
  char* const line = string_get_buffer(tag_line);
  switch (line[0]) {
    case FASTA_TAG_BEGIN:
      sequence->has_qualities = false;
      break;
    case FASTQ_TAG_BEGIN:
      sequence->has_qualities = true;
      break;
    default:
      input_fasta_parser_prompt_error(buffered_input,FASTA_ERROR_TAG_BEGINNING);
      return INPUT_STATUS_FAIL;
      break;
  }
  // Delimit the tag length
  uint64_t i = 1;
  while (i < line_length && line[i]!=TAB && line[i]!=SPACE) ++i; // Read SPACE or TAB
  string_set_buffer(&sequence->tag,line+1,i-1);
  // Chomp sequence-end information
  sequence->end_info = input_fasta_parse_sequence_tag_chomp_end_info(&sequence->tag);
  // Parse extra info
  string_clear(&sequence->extra_tag);
  input_fasta_parse_sequence_tag_extra(sequence,line+i,line_length-i);
  // Return
  return INPUT_STATUS_OK;
}
int input_fasta_parse_sequence_read(
    buffered_input_file_t* const buffered_input,
    sequence_t* const sequence,
    const bool correct_sequence) {
  // Read READ line
  string_t* const read = &sequence->read;
  buffered_input_file_get_line(buffered_input,read);
  // Check/Normalize READ
  const uint64_t line_length = string_get_length(read)-1; // No EOL
  char* const line = string_get_buffer(read);
  bool bad_character_warn = false;
  uint64_t i = 0;
  for (i=0;i<line_length-1;++i) {
    const char character = line[i];
    if (!is_dna(character)) {
      if (!is_iupac_code(character)) {
        if (!correct_sequence) {
          input_fasta_parser_prompt_error(buffered_input,FASTA_ERROR_READ_BAD_CHARACTER);
          return INPUT_STATUS_FAIL;
        }
        if (!bad_character_warn) {
          bad_character_warn = true;
          gem_warn(PARSE_FASTQ_READ_BAD_CHARACTER,
              buffered_input_file_get_file_name(buffered_input),
              buffered_input_file_get_current_line_num(buffered_input)-1,"replaced with 'N'");
        }
      }
      line[i] = dna_normalized(character);
    }
  }
  string_set_length(&sequence->read,line_length);
  return INPUT_STATUS_OK;
}
int input_fasta_parse_sequence_qualities(
    buffered_input_file_t* const buffered_input,
    sequence_t* const sequence) {
  // Read QUALITIES line
  string_t* const qualities = &sequence->qualities;
  buffered_input_file_get_line(buffered_input,qualities);
  // Check/Normalize QUALITIES
  const uint64_t line_length = string_get_length(qualities)-1; // No EOL
  char* const line = string_get_buffer(qualities);
  uint64_t i = 0;
  for (i=0;i<line_length-1;++i) {
    const char character = line[i];
    if (gem_expect_false(!SEQUENCE_QUALITY_IS_VALID(character))) {
      input_fasta_parser_prompt_error(buffered_input,FASTA_ERROR_QUALITIES_BAD_CHARACTER);
      return INPUT_STATUS_FAIL;
    }
  }
  string_set_length(&sequence->qualities,line_length);
  // Return
  return INPUT_STATUS_OK;
}
/*
 * High Level Parsers
 */
int input_fasta_parse_sequence_components(
    buffered_input_file_t* const buffered_input,
    sequence_t* const sequence) {
  // Parse TAG
  int status;
  status = input_fasta_parse_sequence_tag(buffered_input,sequence);
  if (status!=INPUT_STATUS_OK) return status;
  // Parse READ
  status = input_fasta_parse_sequence_read(buffered_input,sequence,true);
  if (status!=INPUT_STATUS_OK) return status;
  // Parse QUALITIES
  if (sequence->has_qualities) {
    // Read line (FASTQ comments)
    buffered_input_file_get_line(buffered_input,&sequence->qualities);
    const char* const buffer = string_get_buffer(&sequence->qualities);
    if (buffer[0]!=FASTQ_SEP) { // Check '+'
      input_fasta_parser_prompt_error(buffered_input,FASTA_ERROR_SEPARATOR_BAD_CHARACTER);
      return INPUT_STATUS_FAIL;
    }
    // Parse qualities string
    status = input_fasta_parse_sequence_qualities(buffered_input,sequence);
    if (status!=INPUT_STATUS_OK) return status;
    // Check lengths
    if (string_get_length(&sequence->read) != string_get_length(&sequence->qualities)) {
      input_fasta_parser_prompt_error(buffered_input,FASTA_ERROR_LENGTHS);
      string_clear(&sequence->qualities);
    }
  }
  return INPUT_STATUS_OK;
}
int input_fasta_parse_sequence(
    buffered_input_file_t* const buffered_input,
    sequence_t* const sequence,
    const bool check_input_buffer) {
  PROFILE_START(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
  // Buffer handling
  int status;
  if (check_input_buffer && buffered_input_file_eob(buffered_input)) {
    status = buffered_input_file_reload(buffered_input,0);
    if (status!=INPUT_STATUS_OK) {
      PROFILE_STOP(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
      return status;
    }
  }
  // Parse FASTQ-record components
  sequence_clear(sequence); // Clear read
  status = input_fasta_parse_sequence_components(buffered_input,sequence);
  if (status!=INPUT_STATUS_OK) {
    PROFILE_STOP(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
    return status;
  }
  PROFILE_STOP(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
  return INPUT_STATUS_OK;
}
/*
 * Display
 */
void input_fasta_parser_prompt_error(
    buffered_input_file_t* const buffered_input,
    const error_code_t error_code) {
  // Display textual error msg
  const char* const file_name = buffered_input_file_get_file_name(buffered_input);
  const uint64_t line_no = buffered_input_file_get_current_line_num(buffered_input)-1;
  switch (error_code) {
    case 0:
      /* No error */
      break;
    case FASTA_ERROR_TAG_BEGINNING:
      gem_fatal_error(PARSE_FASTQ_TAG_BEGINNING,file_name,line_no);
      break;
    case FASTA_ERROR_SEPARATOR_BAD_CHARACTER:
      gem_fatal_error(PARSE_FASTQ_SEPARATOR_BAD_CHARACTER,file_name,line_no);
      break;
    case FASTA_ERROR_QUALITIES_BAD_CHARACTER:
      gem_fatal_error(PARSE_FASTQ_QUALITIES_BAD_CHARACTER,file_name,line_no);
      break;
    case FASTA_ERROR_LENGTHS:
      gem_fatal_error(PARSE_FASTQ_LENGTHS,file_name,line_no);
      break;
    default:
      gem_fatal_error(PARSE_FASTQ,file_name,line_no);
      break;
  }
}
