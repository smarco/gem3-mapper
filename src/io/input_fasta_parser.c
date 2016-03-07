/*
 * PROJECT: GEMMapper
 * FILE: input_fastq_parser.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for FASTA/FASTQ/MULTIFASTA files
 */

#include "io/input_fasta_parser.h"
#include "io/input_parser.h"
#include "data_structures/dna_text.h"

/*
 * Error codes
 */
//#define FASTA_ERROR_WRONG_FILE_FORMAT 10   /**/
#define FASTA_ERROR_PREMATURE_EOF 10
#define FASTA_ERROR_TAG_BEGINNING 20
#define FASTA_ERROR_READ_BAD_CHARACTER 30
#define FASTA_ERROR_SEPARATOR_BAD_CHARACTER 40
#define FASTA_ERROR_QUALITIES_BAD_CHARACTER 50
#define FASTA_ERROR_LENGTHS 60

/*
 * Error Messages
 */
#define GEM_ERROR_PARSE_FASTQ "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64")"
#define GEM_ERROR_PARSE_FASTQ_PREMATURE_EOF "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Premature end-of-file"
#define GEM_ERROR_PARSE_FASTQ_READ_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64"). Input read sequence contains bad character %s"
#define GEM_ERROR_PARSE_FASTQ_QUALITIES_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64"). Input read qualities contains invalid quality value"
#define GEM_ERROR_PARSE_FASTQ_TAG_BEGINNING "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Beginning Symbol ('>' or '@') not found. Bad syntax"
#define GEM_ERROR_PARSE_FASTQ_SEPARATOR_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Wrong separator symbol ('+'). Bad syntax"
#define GEM_ERROR_PARSE_FASTQ_LENGTHS "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Different read and qualities lengths"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Error handler
 */
void input_fasta_parser_prompt_error(
    buffered_input_file_t* const buffered_fasta_input,
    uint64_t line_num,
    uint64_t column_pos,
    const error_code_t error_code) {
  // Display textual error msg
  const char* const file_name = input_file_get_file_name(buffered_fasta_input->input_file);
  switch (error_code) {
    case 0: /* No error */ break;
    case FASTA_ERROR_PREMATURE_EOF: gem_error(PARSE_FASTQ_PREMATURE_EOF,file_name,line_num); break;
    case FASTA_ERROR_TAG_BEGINNING: gem_error(PARSE_FASTQ_TAG_BEGINNING,file_name,line_num); break;
    case FASTA_ERROR_SEPARATOR_BAD_CHARACTER: gem_error(PARSE_FASTQ_SEPARATOR_BAD_CHARACTER,file_name,line_num); break;
    case FASTA_ERROR_LENGTHS: gem_error(PARSE_FASTQ_LENGTHS,file_name,line_num); break;
    default: gem_error(PARSE_FASTQ,file_name,line_num,column_pos); break;
  }
}
/*
 * FASTQ/FASTA file. Skip record
 *   - Here, the trick is to try to re-synch with FASTQ/FASTA records after a syntax error
 */
#define FASTA_SEEMS_TAG_COND(input_buffer) \
  (input_buffer_eob(input_buffer) || \
   input_buffer->cursor[0]==FASTQ_TAG_BEGIN || \
   input_buffer->cursor[0]==FASTA_TAG_BEGIN)
#define FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(input_buffer) \
  while (!input_buffer_eob(input_buffer) && \
         input_buffer->cursor[0]!=FASTQ_TAG_BEGIN && \
         input_buffer->cursor[0]!=FASTA_TAG_BEGIN) { \
    input_buffer_skip_line(input_buffer); \
  }
void input_fasta_parser_next_record(
    buffered_input_file_t* const buffered_fasta_input,
    char* const line_start) {
  input_buffer_t* const input_buffer = buffered_fasta_input->input_buffer;
  // Back to the beginning
  input_buffer->cursor = line_start;
  // Skip first line(s)
  if (input_buffer_eob(input_buffer)) return; // oops!!
  if (input_buffer->cursor[0]==FASTA_TAG_BEGIN) { // Seems FASTA
    input_buffer_skip_line(input_buffer); // Skip TAG
    if (FASTA_SEEMS_TAG_COND(input_buffer)) return; // Let's hope we're done
    input_buffer_skip_line(input_buffer); // Skip READ
    if (FASTA_SEEMS_TAG_COND(input_buffer)) return; // Let's hope we're done
    FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(input_buffer);
  } else if (input_buffer->cursor[0]==FASTQ_TAG_BEGIN) { // Seems FASTQ
    input_buffer_skip_line(input_buffer); // Skip TAG
    if (FASTA_SEEMS_TAG_COND(input_buffer)) return; // Let's hope we're done
    input_buffer_skip_line(input_buffer); // Skip READ
    if (FASTA_SEEMS_TAG_COND(input_buffer)) return; // Let's hope we're done
    input_buffer_skip_line(input_buffer); // Skip '+'
    if (FASTA_SEEMS_TAG_COND(input_buffer)) return; // Let's hope we're done
    input_buffer_skip_line(input_buffer); // Skip QUALS
    if (FASTA_SEEMS_TAG_COND(input_buffer)) return; // Let's hope we're done
    FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(input_buffer);
  } else { // Weird
    FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(input_buffer);
  }
  return; // Let's hope we can get to synchronize at some point
}
/*
 * Accessors
 */
bool input_fasta_is_fasta(input_file_t* const input_file) {
  return input_file->file_format==FASTA && input_file->fasta==F_FASTA;
}
bool input_fasta_is_fastq(input_file_t* const input_file) {
  return input_file->file_format==FASTA && input_file->fasta==F_FASTQ;
}
/*
 * FASTQ/FASTA format. Basic building block for parsing
 */
error_code_t ifp_parse_tag(
    input_buffer_t* const input_buffer,
    string_t* const tag,
    sequence_attributes_t* const attributes,
    bool* const has_qualities) {
  // Beginning character
  char** const text_line = input_buffer_get_cursor(input_buffer);
  switch (**text_line) {
    case FASTA_TAG_BEGIN:
      *has_qualities=false;
      break;
    case FASTQ_TAG_BEGIN:
      *has_qualities=true;
      break;
    default:
      return FASTA_ERROR_TAG_BEGINNING;
      break;
  }
  PARSER_NEXT_CHAR(text_line);
  // Parse TAG
  input_text_parse_tag(text_line,tag,attributes);
  if (PARSER_IS_EOL(text_line)) PARSER_NEXT_CHAR(text_line);
  return 0;
}
error_code_t ifp_parse_read(
    input_file_t* const input_file,
    input_buffer_t* const input_buffer,
    string_t* const read,
    const bool strictly_normalized,
    const bool correct_sequence) {
  if (input_buffer_eob(input_buffer)) return FASTA_ERROR_PREMATURE_EOF;
  char** const text_line = input_buffer_get_cursor(input_buffer);
  char* const read_begin = *text_line;
  bool bad_character_warn = false;
  while (!PARSER_IS_EOL(text_line)) {
    const char character = **text_line;
    if (gem_expect_false(!is_dna(character))) {
      if (!is_iupac_code(**text_line)) {
        if (!correct_sequence) return FASTA_ERROR_READ_BAD_CHARACTER;
        if (!bad_character_warn) {
          bad_character_warn = true;
          gem_warn(PARSE_FASTQ_READ_BAD_CHARACTER,
              input_file_get_file_name(input_file),input_buffer->current_line_num,
              (uint64_t)(*text_line-read_begin),"replaced with 'N'");
        }
      }
      **text_line = (strictly_normalized) ? dna_strictly_normalized(**text_line) : dna_normalized(**text_line);
    }
    PARSER_NEXT_CHAR(text_line);
  }
  // Set the read
  const uint64_t length = (*text_line-read_begin);
  string_set_buffer(read,read_begin,length);
  PARSER_NEXT_CHAR(text_line);
  return 0;
}
error_code_t ifp_parse_qualities(
    input_file_t* const input_file,
    input_buffer_t* const input_buffer,
    string_t* const qualities,
    const bool try_recovery) {
  if (input_buffer_eob(input_buffer)) return FASTA_ERROR_PREMATURE_EOF;
  char** const text_line = input_buffer_get_cursor(input_buffer);
  char* const quals_begin = *text_line;
  while (!PARSER_IS_EOL(text_line)) {
    if (gem_expect_false(!SEQUENCE_QUALITY_IS_VALID(**text_line))) {
      if (try_recovery) {
        gem_warn(PARSE_FASTQ_QUALITIES_BAD_CHARACTER,input_file_get_file_name(input_file),
            (uint64_t)input_buffer->current_line_num,(uint64_t)(*text_line-quals_begin));
        string_clear(qualities);
        return 0;
      }
      return FASTA_ERROR_QUALITIES_BAD_CHARACTER;
    }
    PARSER_NEXT_CHAR(text_line);
  }
  // Set the qualities
  const uint64_t length = (*text_line-quals_begin);
  string_set_buffer(qualities,quals_begin,length);
  PARSER_NEXT_CHAR(text_line);
  return 0;
}
/*
 * High Level Parsers
 */
error_code_t ifp_parse_sequence(
    input_file_t* const input_file,
    input_buffer_t* const input_buffer,
    sequence_t* const seq_read,
    const bool strictly_normalized,
    const bool try_recovery) {
  // Parse TAG
  error_code_t error_code;
  bool has_qualities;
  if ((error_code=ifp_parse_tag(input_buffer,&seq_read->tag,&seq_read->attributes,&has_qualities))) return error_code;
  // Parse READ
  if ((error_code=ifp_parse_read(input_file,input_buffer,&seq_read->read,true,true))) return error_code;
  // Parse QUALITIES
  if (has_qualities) {
    // Skip '+'
    char** const text_line = input_buffer_get_cursor(input_buffer);
    if (**text_line!=FASTQ_SEP) {
      return input_buffer_eob(input_buffer) ? FASTA_ERROR_PREMATURE_EOF : FASTA_ERROR_SEPARATOR_BAD_CHARACTER;
    }
    PARSER_SKIP_LINE(text_line); PARSER_NEXT_CHAR(text_line);
    // Parse qualities string
    if ((error_code=ifp_parse_qualities(input_file,input_buffer,&seq_read->qualities,try_recovery))) return error_code;
    // Check lengths
    if (gem_expect_false(string_get_length(&seq_read->read) != string_get_length(&seq_read->qualities))) {
      if (!try_recovery) return FASTA_ERROR_LENGTHS;
      string_clear(&seq_read->qualities);
    }
    input_buffer->current_line_num += 4;
  } else {
    input_buffer->current_line_num += 2;
  }
  return 0;
}
error_code_t input_fasta_parse_sequence(
    buffered_input_file_t* const buffered_fasta_input,
    sequence_t* const seq_read,
    const bool strictly_normalized,
    const bool try_recovery,
    const bool check_input_buffer) {
  PROFILE_START(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
  error_code_t error_code;
  if (check_input_buffer && buffered_input_file_eob(buffered_fasta_input)) {
    if ((error_code=buffered_input_file_reload__dump_attached(buffered_fasta_input,0))!=INPUT_STATUS_OK) {
      PROFILE_STOP(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
      return error_code;
    }
  }
  // Parse read
  input_buffer_t* const input_buffer = buffered_fasta_input->input_buffer;
  char* const line_start = input_buffer->cursor;
  const uint64_t line_num = input_buffer->current_line_num;
  sequence_clear(seq_read); // Prepare read
  if ((error_code=ifp_parse_sequence(
      buffered_fasta_input->input_file,input_buffer,seq_read,strictly_normalized,try_recovery))) {
    const uint64_t column_pos = input_buffer->cursor-line_start;
    input_fasta_parser_prompt_error(buffered_fasta_input,line_num,column_pos,error_code);
    input_fasta_parser_next_record(buffered_fasta_input,line_start);
    PROFILE_STOP(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
    return INPUT_STATUS_FAIL;
  }
  PROFILE_STOP(GP_INPUT_FASTA_PARSE_SEQUENCE,PROFILE_LEVEL);
  return INPUT_STATUS_OK;
}
