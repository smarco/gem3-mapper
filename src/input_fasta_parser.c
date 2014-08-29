/*
 * PROJECT: GEMMapper
 * FILE: input_fastq_parser.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for FASTA/FASTQ/MULTIFASTA files
 */

#include "input_fasta_parser.h"
#include "input_parser.h"
#include "dna_string.h"

// Constants
#define FASTA_BUFFER_NUM_LINES (2*4*NUM_LINES_5K) /* 2l-Paired x 4l-FASTQRecord x 5K-BufferSize */

/*
 * FASTQ/FASTA File Format test
 */
/*
#define FASTA_TEST_FASTA_SKIP_LINE() \
  while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL && buffer[buffer_pos]!=DOS_EOL) ++buffer_pos; \
  INPUT_FILE_SKIP_EOL(buffer,buffer_pos)
GEM_INLINE bool input_fasta_parser_test_fastq(
    char* const file_name,const uint64_t line_num,char* const buffer,const uint64_t buffer_size,
    fasta_file_format_t* const fasta_file_format,const bool show_errors) {
  // Check trivial case
  if (gem_expect_false(buffer_size==0)) return true;
  // Test TAG {fastq,fasta}
  bool expect_qualities;
  switch (buffer[0]) {
    case FASTA_TAG_BEGIN:
      expect_qualities=false;
      break;
    case FASTQ_TAG_BEGIN:
      expect_qualities=true;
      break;
    default:
      return false;
      break;
  }
  uint64_t buffer_pos=1;
  FASTA_TEST_FASTA_SKIP_LINE(); // Skip the rest of the TAG
  if (buffer_pos==buffer_size) return false;
  // Test READ
  const uint64_t read_start_pos = buffer_pos;
  while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL && buffer[buffer_pos]!=DOS_EOL) {
    if (gem_expect_false(!is_iupac_code(buffer[buffer_pos]))) return false;
    ++buffer_pos;
  }
  INPUT_FILE_SKIP_EOL(buffer,buffer_pos);
  if (buffer_pos==buffer_size) return expect_qualities;
  const uint64_t read_length = buffer_pos-read_start_pos;
  // Test SEPARATOR
  switch (buffer[buffer_pos]) {
    case FASTQ_SEP: // '+'
      if (!expect_qualities) return false;
      // Skip remarks
      ++buffer_pos;
      FASTA_TEST_FASTA_SKIP_LINE();
      if (buffer_pos==buffer_size) return false;
      // Found regular FASTQ !
      *fasta_file_format = F_FASTQ;
      return true;
      break;
    case FASTA_TAG_BEGIN: // '>'
      if (expect_qualities) return false;
      // Found regular FASTA !
      *fasta_file_format = F_FASTA;
      return true;
      break;
    default:
      return false;
  }
}
GEM_INLINE bool input_file_test_fasta(
    input_file_t* const input_file,fasta_file_format_t* const fasta_file_format,const bool show_errors) {
  INPUT_FILE_CHECK(input_file);
  return (input_fasta_parser_test_fastq(input_file_get_file_name(input_file),input_file->processed_lines+1,
      (char*)input_file->file_buffer,input_file->buffer_size,fasta_file_format,show_errors));
}
GEM_INLINE error_code_t input_fasta_parser_check_fastq_file_format(buffered_input_file_t* const buffered_fasta_input) {
  input_file_t* const input_file = buffered_fasta_input->input_file;
  if (expect_false(input_file->file_format==FILE_FORMAT_UNKNOWN)) { // Unknown
    fasta_file_format_t fasta_file_format;
    if (!input_fasta_parser_test_fastq(
        input_file_get_file_name(input_file),buffered_fasta_input->current_line_num,
        vector_get_mem(buffered_fasta_input->block_buffer,char),
        vector_get_used(buffered_fasta_input->block_buffer),&fasta_file_format,true)) {
      return IFP_PE_WRONG_FILE_FORMAT;
    }
    input_file->file_format = FASTA;
    input_file->fasta = fasta_file_format;
  } else if (expect_false(input_file->file_format!=FASTA)) {
    return IFP_PE_WRONG_FILE_FORMAT;
  }
  return 0;
}
*/


/*
 * Error handler
 */
GEM_INLINE void input_fasta_parser_prompt_error(
    buffered_input_file_t* const buffered_fasta_input,
    uint64_t line_num,uint64_t column_pos,const error_code_t error_code) {
  // Display textual error msg
  const char* const file_name = input_file_get_file_name(buffered_fasta_input->input_file);
  switch (error_code) {
    case 0: /* No error */ break;
    case FASTA_ERROR_TAG_BEGINNING: gem_error(FASTA_TAG_BEGINNING,file_name,line_num); break;
    case FASTA_ERROR_SEPARATOR_BAD_CHARACTER: gem_error(FASTA_SEPARATOR_BAD_CHARACTER,file_name,line_num); break;
    case FASTA_ERROR_LENGTHS: gem_error(FASTA_LENGTHS,file_name,line_num); break;
    default:
      gem_error(PARSE_FASTA,file_name,line_num,column_pos);
      break;
  }
}
/*
 * FASTQ/FASTA file. Skip record
 *   - Here, the trick is to try to re-synch with FASTQ/FASTA records after a syntax error
 */
#define FASTA_SEEMS_TAG_COND(buffered_fasta_input) \
  (buffered_input_file_eob(buffered_fasta_input) || \
   buffered_fasta_input->cursor[0]==FASTQ_TAG_BEGIN || \
   buffered_fasta_input->cursor[0]==FASTA_TAG_BEGIN)
#define FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(buffered_fasta_input) \
  while (!buffered_input_file_eob(buffered_fasta_input) && \
      buffered_fasta_input->cursor[0]!=FASTQ_TAG_BEGIN && \
      buffered_fasta_input->cursor[0]!=FASTA_TAG_BEGIN) { \
    buffered_input_file_skip_line(buffered_fasta_input); \
  }
GEM_INLINE void input_fasta_parser_next_record(buffered_input_file_t* const buffered_fasta_input,char* const line_start) {
  BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  GEM_CHECK_NULL(line_start);
  // Back to the beginning
  buffered_fasta_input->cursor = line_start;
  // Skip first line(s)
  if (buffered_input_file_eob(buffered_fasta_input)) return; // oops!!
  if (buffered_fasta_input->cursor[0]==FASTA_TAG_BEGIN) { // Seems FASTA
    buffered_input_file_skip_line(buffered_fasta_input); // Skip TAG
    if (FASTA_SEEMS_TAG_COND(buffered_fasta_input)) return; // Let's hope we're done
    buffered_input_file_skip_line(buffered_fasta_input); // Skip READ
    if (FASTA_SEEMS_TAG_COND(buffered_fasta_input)) return; // Let's hope we're done
    FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(buffered_fasta_input);
  } else if (buffered_fasta_input->cursor[0]==FASTQ_TAG_BEGIN) { // Seems FASTQ
    buffered_input_file_skip_line(buffered_fasta_input); // Skip TAG
    if (FASTA_SEEMS_TAG_COND(buffered_fasta_input)) return; // Let's hope we're done
    buffered_input_file_skip_line(buffered_fasta_input); // Skip READ
    if (FASTA_SEEMS_TAG_COND(buffered_fasta_input)) return; // Let's hope we're done
    buffered_input_file_skip_line(buffered_fasta_input); // Skip '+'
    if (FASTA_SEEMS_TAG_COND(buffered_fasta_input)) return; // Let's hope we're done
    buffered_input_file_skip_line(buffered_fasta_input); // Skip QUALS
    if (FASTA_SEEMS_TAG_COND(buffered_fasta_input)) return; // Let's hope we're done
    FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(buffered_fasta_input);
  } else { // Weird
    FASTA_SKIP_LINES_UNTIL_SYNCH_TAG(buffered_fasta_input);
  }
  return; // Let's hope we can get to synchronize at some point
}
/*
 * Accessors
 */
GEM_INLINE bool input_fasta_is_fasta(input_file_t* const input_file) {
  return input_file->file_format==FASTA && input_file->fasta==F_FASTA;
}
GEM_INLINE bool input_fasta_is_fastq(input_file_t* const input_file) {
  return input_file->file_format==FASTA && input_file->fasta==F_FASTQ;
}
/*
 * FASTQ/FASTA format. Basic building block for parsing
 */
GEM_INLINE error_code_t ifp_parse_tag(char** const text_line,string_t* const tag,sequence_attributes_t* const attributes,bool* const has_qualities) {
  TEXT_LINE_CHECK(text_line);
  STRING_CHECK(tag);
  // Check begin
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
GEM_INLINE error_code_t ifp_parse_read(
    buffered_input_file_t* const buffered_fasta_input,string_t* const read,
    const bool strictly_normalized,const bool correct_sequence) {
  char** const text_line = buffered_input_file_get_text_line(buffered_fasta_input);
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
              input_file_get_file_name(buffered_fasta_input->input_file),
              buffered_fasta_input->current_line_num,(uint64_t)(*text_line-read_begin),"replaced with 'N'");
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
GEM_INLINE error_code_t ifp_parse_qualities(
    buffered_input_file_t* const buffered_fasta_input,string_t* const qualities,const bool try_recovery) {
  char** const text_line = buffered_input_file_get_text_line(buffered_fasta_input);
  char* const quals_begin = *text_line;
  while (!PARSER_IS_EOL(text_line)) {
    if (gem_expect_false(!SEQUENCE_QUALITY_IS_VALID(**text_line))) {
      if (try_recovery) {
        gem_warn(PARSE_FASTQ_QUALITIES_BAD_CHARACTER,
            input_file_get_file_name(buffered_fasta_input->input_file),
            buffered_fasta_input->current_line_num,*text_line-quals_begin);
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
GEM_INLINE error_code_t ifp_parse_sequence(
    buffered_input_file_t* const buffered_fasta_input,sequence_t* const seq_read,
    const bool strictly_normalized,const bool try_recovery) {
  // Parse TAG
  error_code_t error_code;
  bool has_qualities;
  if ((error_code=ifp_parse_tag(buffered_input_file_get_text_line(buffered_fasta_input),
      seq_read->tag,&(seq_read->attributes),&has_qualities))) return error_code;
  // Parse READ
  if ((error_code=ifp_parse_read(buffered_fasta_input,seq_read->read,true,true))) return error_code;
  // Parse QUALITIES
  if (has_qualities) {
    // Skip '+'
    char** const text_line = buffered_input_file_get_text_line(buffered_fasta_input);
    if (**text_line!=FASTQ_SEP) return FASTA_ERROR_SEPARATOR_BAD_CHARACTER;
    PARSER_SKIP_LINE(text_line); PARSER_NEXT_CHAR(text_line);
    // Parse qualities string
    if ((error_code=ifp_parse_qualities(buffered_fasta_input,seq_read->qualities,try_recovery))) return error_code;
    // Check lengths
    if (gem_expect_false(string_get_length(seq_read->read) != string_get_length(seq_read->qualities))) {
      if (!try_recovery) return FASTA_ERROR_LENGTHS;
      string_clear(seq_read->qualities);
    }
    buffered_fasta_input->current_line_num += 4;
  } else {
    buffered_fasta_input->current_line_num += 2;
  }
  return 0;
}
GEM_INLINE error_code_t input_fasta_parse_sequence_(
    buffered_input_file_t* const buffered_fasta_input,sequence_t* const seq_read,
    const bool strictly_normalized,const bool try_recovery) {
  BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  SEQUENCE_CHECK(seq_read);
  // Prepare read
  sequence_clear(seq_read);
  // Parse read
  error_code_t error_code;
  char* const line_start = buffered_fasta_input->cursor;
  const uint64_t line_num = buffered_fasta_input->current_line_num;
  if ((error_code=ifp_parse_sequence(buffered_fasta_input,seq_read,strictly_normalized,try_recovery))) {
    const uint64_t column_pos = buffered_fasta_input->cursor-line_start;
    input_fasta_parser_prompt_error(buffered_fasta_input,line_num,column_pos,error_code);
    input_fasta_parser_next_record(buffered_fasta_input,line_start);
    return INPUT_STATUS_FAIL;
  }
  return INPUT_STATUS_OK;
}
GEM_INLINE error_code_t input_fasta_parse_sequence(
    buffered_input_file_t* const buffered_fasta_input,sequence_t* const seq_read,
    const bool strictly_normalized,const bool try_recovery) {
  BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  // Check input file
  error_code_t error_code;
  if (buffered_input_file_eob(buffered_fasta_input)) { // Check the end_of_block. Reload buffer if needed
    if ((error_code=buffered_input_file_reload(buffered_fasta_input,FASTA_BUFFER_NUM_LINES))!=INPUT_STATUS_OK) return error_code;
  }
  //  // Check file format
  //  if (input_fasta_parser_check_fastq_file_format(buffered_fasta_input)) {
  //    gem_fatal_error(PARSE_FASTQ_BAD_FILE_FORMAT,
  //        input_file_get_file_name(input_file),
  //        buffered_input_file_get_current_line(buffered_fasta_input));
  //    return INPUT_STATUS_FAIL;
  //  }
  // Parse read
  return input_fasta_parse_sequence_(buffered_fasta_input,seq_read,strictly_normalized,try_recovery);
}
