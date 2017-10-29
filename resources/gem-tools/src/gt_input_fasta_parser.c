/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_fastq_parser.c
 * DATE: 17/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_input_fasta_parser.h"

// Constants
#define GT_IFP_NUM_LINES (2/*Paired*/*4/*4_lines_per_record*/*GT_NUM_LINES_5K)
#define GT_IFP_MULTIFASTA_NUM_LINES GT_NUM_LINES_2M

#define GT_IFP_FASTA_TAG_BEGIN '>'
#define GT_IFP_FASTQ_TAG_BEGIN '@'
#define GT_IFP_FASTQ_SEP '+'

/*
 * FASTQ/FASTA File Format test
 */
#define GT_IFP_TEST_FASTA_SKIP_LINE() \
  while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL && buffer[buffer_pos]!=DOS_EOL) ++buffer_pos; \
  if (buffer[buffer_pos]==DOS_EOL) ++buffer_pos; \
  if (buffer_pos==buffer_size) return false
GT_INLINE bool gt_input_fasta_parser_test_fastq(
    char* const file_name,const uint64_t line_num,char* const buffer,const uint64_t buffer_size,
    gt_fasta_file_format* const fasta_file_format,const bool show_errors) {
  // Get TAG {fastq,fasta}
  bool expect_qualities;
  switch (buffer[0]) {
    case GT_IFP_FASTA_TAG_BEGIN:
      expect_qualities=false;
      break;
    case GT_IFP_FASTQ_TAG_BEGIN:
      expect_qualities=true;
      break;
    default:
      return false;
      break;
  }
  // Skip TAG
  uint64_t buffer_pos=1;
  GT_IFP_TEST_FASTA_SKIP_LINE();
  ++buffer_pos; // Skip EOL
  // Skip Read
  const uint64_t read_start_pos = buffer_pos;
  while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL && buffer[buffer_pos]!=DOS_EOL) {
    if (gt_expect_false(!gt_is_iupac_code(buffer[buffer_pos]))) return false;
    ++buffer_pos;
  }
  if (buffer[buffer_pos]==DOS_EOL) ++buffer_pos;
  if (buffer_pos==buffer_size) return expect_qualities;
  const uint64_t read_length = buffer_pos-read_start_pos;
  ++buffer_pos;
  // Read next character
  switch (buffer[buffer_pos]) {
    case GT_IFP_FASTQ_SEP: // '+'
      if (!expect_qualities) return false;
      ++buffer_pos;
      // Skip remarks
      GT_IFP_TEST_FASTA_SKIP_LINE();
      ++buffer_pos;
      // Read qualities
      const uint64_t qual_start_pos = buffer_pos;
      while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL && buffer[buffer_pos]!=DOS_EOL) ++buffer_pos;
      if (buffer[buffer_pos]==DOS_EOL) ++buffer_pos;
      const uint64_t qual_length = buffer_pos-qual_start_pos;
      // Check lengths
      if (qual_length!=read_length) return false;
      // Found regular FASTQ !
      fasta_file_format->fasta_format = F_FASTQ;
      return true;
      break;
    case GT_IFP_FASTA_TAG_BEGIN: // '>'
      if (expect_qualities) return false;
      // Found regular FASTA !
      fasta_file_format->fasta_format = F_FASTA;
      return true;
      break;
    default:
      if (!expect_qualities && gt_is_iupac_code(buffer[buffer_pos])) {
        // Found regular MULTI_FASTA !
        fasta_file_format->fasta_format = F_MULTI_FASTA;
        return true;
      }
      return false;
  }
}
GT_INLINE bool gt_input_file_test_fasta(
    gt_input_file* const input_file,gt_fasta_file_format* const fasta_file_format,const bool show_errors) {
  GT_INPUT_FILE_CHECK(input_file);
  return (gt_input_fasta_parser_test_fastq(
      gt_input_file_get_file_name(input_file),input_file->processed_lines+1,
      (char*)input_file->file_buffer,input_file->buffer_size,fasta_file_format,show_errors));
}
GT_INLINE gt_status gt_input_fasta_parser_check_fastq_file_format(gt_buffered_input_file* const buffered_fasta_input) {
  gt_input_file* const input_file = buffered_fasta_input->input_file;
  if (gt_expect_false(input_file->file_format==FILE_FORMAT_UNKNOWN)) { // Unknown
    gt_fasta_file_format fasta_file_format;
    if (!gt_input_fasta_parser_test_fastq(
        gt_input_file_get_file_name(input_file),buffered_fasta_input->current_line_num,
        gt_vector_get_mem(buffered_fasta_input->block_buffer,char),
        gt_vector_get_used(buffered_fasta_input->block_buffer),&fasta_file_format,true)) {
      return GT_IFP_PE_WRONG_FILE_FORMAT;
    }
    input_file->file_format = FASTA;
    input_file->fasta_type = fasta_file_format;
  } else if (gt_expect_false(input_file->file_format!=FASTA)) {
    return GT_IFP_PE_WRONG_FILE_FORMAT;
  }
  return 0;
}

/*
 * FASTQ/FASTA File basics
 */
/* Error handler */
GT_INLINE void gt_input_fasta_parser_prompt_error(
    gt_buffered_input_file* const buffered_fasta_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code) {
  // Display textual error msg
  const char* const file_name = gt_input_file_get_file_name(buffered_fasta_input->input_file);
  switch (error_code) {
    case 0: /* No error */ break; // TODO
    // case GT_IMP_PE_WRONG_FILE_FORMAT: gt_error(PARSE_MAP_BAD_FILE_FORMAT,file_name,line_num); break;
    default:
      gt_error(PARSE_FASTA,file_name,line_num,column_pos);
      break;
  }
}
#define GT_IFP_SEEMS_TAG_COND() \
  (gt_buffered_input_file_eob(buffered_fasta_input) || \
   buffered_fasta_input->cursor[0]==GT_IFP_FASTQ_TAG_BEGIN || \
   buffered_fasta_input->cursor[0]==GT_IFP_FASTA_TAG_BEGIN)
#define GT_IFP_SKIP_LINES_UNTIL_SYNCH_TAG() \
  while (!gt_buffered_input_file_eob(buffered_fasta_input) && \
      buffered_fasta_input->cursor[0]!=GT_IFP_FASTQ_TAG_BEGIN && \
      buffered_fasta_input->cursor[0]!=GT_IFP_FASTA_TAG_BEGIN) { \
    gt_buffered_input_file_skip_line(buffered_fasta_input); \
  }
/* FASTQ/FASTA file. Skip record */
GT_INLINE void gt_input_fasta_parser_next_record(gt_buffered_input_file* const buffered_fasta_input,char* const line_start) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  GT_NULL_CHECK(line_start);
  /*
   * Here, the trick is to try to resynch with FASTQ/FASTA records after a syntax error
   */
  // Back to the beginning
  buffered_fasta_input->cursor = line_start;
  // Skip first line(s)
  if (gt_buffered_input_file_eob(buffered_fasta_input)) return; // oops!!
  if (buffered_fasta_input->cursor[0]==GT_IFP_FASTA_TAG_BEGIN) { // Seems FASTA
    gt_buffered_input_file_skip_line(buffered_fasta_input); // Skip TAG
    if (GT_IFP_SEEMS_TAG_COND()) return; // Let's hope we're done
    gt_buffered_input_file_skip_line(buffered_fasta_input); // Skip READ
    if (GT_IFP_SEEMS_TAG_COND()) return; // Let's hope we're done
    GT_IFP_SKIP_LINES_UNTIL_SYNCH_TAG();
  } else if (buffered_fasta_input->cursor[0]==GT_IFP_FASTQ_TAG_BEGIN) { // Seems FASTQ
    gt_buffered_input_file_skip_line(buffered_fasta_input); // Skip TAG
    if (GT_IFP_SEEMS_TAG_COND()) return; // Let's hope we're done
    gt_buffered_input_file_skip_line(buffered_fasta_input); // Skip READ
    if (GT_IFP_SEEMS_TAG_COND()) return; // Let's hope we're done
    gt_buffered_input_file_skip_line(buffered_fasta_input); // Skip '+'
    if (GT_IFP_SEEMS_TAG_COND()) return; // Let's hope we're done
    gt_buffered_input_file_skip_line(buffered_fasta_input); // Skip QUALS
    if (GT_IFP_SEEMS_TAG_COND()) return; // Let's hope we're done
    GT_IFP_SKIP_LINES_UNTIL_SYNCH_TAG();
  } else { // Weird
    GT_IFP_SKIP_LINES_UNTIL_SYNCH_TAG();
  }
  return; // Let's hope we can get to synchronize at some point
}
/* FASTQ/FASTA file. Reload internal buffer */
GT_INLINE gt_status gt_input_fasta_parser_reload_buffer(gt_buffered_input_file* const buffered_fasta_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  // Dump buffer if BOF it attached to input, and get new out block (always FIRST)
  gt_buffered_input_file_dump_attached_buffers(buffered_fasta_input->attached_buffered_output_file);
  // Read new input block
  const uint64_t read_lines =
      gt_buffered_input_file_get_block(buffered_fasta_input,GT_IFP_NUM_LINES);
  if (gt_expect_false(read_lines==0)) return GT_IFP_EOF;
  // Assign block ID
  gt_buffered_input_file_set_id_attached_buffers(buffered_fasta_input->attached_buffered_output_file,buffered_fasta_input->block_id);
  return GT_IFP_OK;
}
/*
 * FASTQ/FASTA format. Basic building block for parsing
 */
GT_INLINE gt_status gt_input_fasta_parse_tag(const char** const text_line,gt_string* const tag,gt_attributes* const attributes) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_STRING_CHECK(tag);
  // Check begin
  if (gt_expect_false(!(**text_line==GT_IFP_FASTA_TAG_BEGIN||**text_line==GT_IFP_FASTQ_TAG_BEGIN))) {
    return GT_IFP_PE_TAG_BAD_BEGINNING;
  }
  GT_NEXT_CHAR(text_line);
  gt_input_parse_tag(text_line,tag,attributes);
  if (GT_IS_EOL(text_line)) GT_NEXT_CHAR(text_line);
  return 0;
}

#define GT_INPUT_FASTQ_PARSE_READ_CHARS(read_begin,length) \
  const char* const read_begin = *text_line; \
  while (!GT_IS_EOL(text_line)) { \
    if (gt_expect_false(!gt_is_iupac_code(**text_line))) return GT_IFP_PE_READ_BAD_CHARACTER; \
    GT_NEXT_CHAR(text_line); \
  } \
  const uint64_t length = (*text_line-read_begin)
GT_INLINE gt_status gt_input_fasta_parse_read_s(const char** const text_line,gt_string* const read) {
  GT_INPUT_FASTQ_PARSE_READ_CHARS(read_begin,length);
  gt_string_append_string(read,read_begin,length);
  GT_NEXT_CHAR(text_line);
  return 0;
}
GT_INLINE gt_status gt_input_fasta_parse_read_sq(const char** const text_line,gt_segmented_sequence* const segmented_sequence) {
  GT_INPUT_FASTQ_PARSE_READ_CHARS(read_begin,length);
  gt_segmented_sequence_append_string(segmented_sequence,read_begin,length);
  GT_NEXT_CHAR(text_line);
  return 0;
}

#define GT_INPUT_FASTA_PARSE_QUALITIES_CHARS(quals_begin,length) \
  const char* const quals_begin = *text_line; \
  while (!GT_IS_EOL(text_line)) { \
    if (gt_expect_false(!gt_is_valid_quality(**text_line))) return GT_IFP_PE_QUALS_BAD_CHARACTER; \
    GT_NEXT_CHAR(text_line); \
  } \
  const uint64_t length = (*text_line-quals_begin)
GT_INLINE gt_status gt_input_fasta_parse_qualities_s(const char** const text_line,gt_string* const read) {
  GT_INPUT_FASTA_PARSE_QUALITIES_CHARS(quals_begin,length);
  gt_string_append_string(read,quals_begin,length);
  GT_NEXT_CHAR(text_line);
  return 0;
}
GT_INLINE gt_status gt_input_fasta_parse_qualities_sq(const char** const text_line,gt_segmented_sequence* const segmented_sequence) {
  GT_INPUT_FASTA_PARSE_QUALITIES_CHARS(quals_begin,length);
  gt_segmented_sequence_append_string(segmented_sequence,quals_begin,length);
  GT_NEXT_CHAR(text_line);
  return 0;
}


/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_ifp_parse_read(
    const char** const text_line,gt_string* const tag,gt_string* const read,
    const bool has_qualitites,gt_string* const qualities,gt_attributes* const attributes) {
  // Parse TAG
  gt_status error_code;
  if ((error_code=gt_input_fasta_parse_tag(text_line,tag,attributes))) return error_code;
  // Parse READ
  if ((error_code=gt_input_fasta_parse_read_s(text_line,read))) return error_code;
  // Parse QUALITIES
  if (has_qualitites) {
    // Skip '+'
    if (**text_line!=GT_IFP_FASTQ_SEP) return GT_IFP_PE_SEPARATOR_BAD_CHARACTER;
    GT_SKIP_LINE(text_line); GT_NEXT_CHAR(text_line);
    // Parse qualities string
    if ((error_code=gt_input_fasta_parse_qualities_s(text_line,qualities))) return error_code;
    // Check lengths
    if (gt_expect_false(gt_string_get_length(read)!=gt_string_get_length(qualities))) return GT_IFP_PE_QUALS_BAD_LENGTH;
  }
  return 0;
}

GT_INLINE gt_status gt_ifp_parse_fasta_fastq_read(
    gt_buffered_input_file* const buffered_fasta_input,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_attributes* const attributes) {
  /*
   * Check input file
   */
  gt_input_file* const input_file = buffered_fasta_input->input_file;
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_fasta_input)) {
    if ((error_code=gt_input_fasta_parser_reload_buffer(buffered_fasta_input))!=GT_IFP_OK) return error_code;
  }
  // Check file format
  if (gt_input_fasta_parser_check_fastq_file_format(buffered_fasta_input)) {
    gt_fatal_error(PARSE_FASTA_BAD_FILE_FORMAT,gt_input_file_get_file_name(input_file),buffered_fasta_input->current_line_num);
    return GT_IFP_FAIL;
  }
  const gt_file_fasta_format fasta_format = input_file->fasta_type.fasta_format;
  if (fasta_format!=F_FASTA && fasta_format!=F_FASTQ) return GT_IFP_PE_WRONG_FILE_FORMAT;
  /*
   * Parse read
   */
  char* const line_start = buffered_fasta_input->cursor;
  const uint64_t line_num = buffered_fasta_input->current_line_num;
  // Parse read

  if ((error_code=gt_ifp_parse_read((const char** const)&buffered_fasta_input->cursor,
      tag,read,fasta_format==F_FASTQ,qualities,attributes))) {
    gt_input_fasta_parser_prompt_error(buffered_fasta_input,line_num,buffered_fasta_input->cursor-line_start,error_code);
    gt_input_fasta_parser_next_record(buffered_fasta_input,line_start);
    return GT_IFP_FAIL;
  }
  // Update total number lines read
  if (fasta_format==F_FASTA) buffered_fasta_input->current_line_num+=2;
  if (fasta_format==F_FASTQ) buffered_fasta_input->current_line_num+=4;
  return GT_IFP_OK;
}


GT_INLINE gt_status gt_input_fasta_parser_get_read(
    gt_buffered_input_file* const buffered_fasta_input,gt_dna_read* const dna_read) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  GT_DNA_READ_CHECK(dna_read);
  // Prepare read
  gt_dna_read_clear(dna_read);
  // Parse FASTA/FASTQ record
  return gt_ifp_parse_fasta_fastq_read(buffered_fasta_input,
      dna_read->tag,dna_read->read,dna_read->qualities,dna_read->attributes);
}
GT_INLINE gt_status gt_input_fasta_parser_get_alignment(
    gt_buffered_input_file* const buffered_fasta_input,gt_alignment* const alignment) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  GT_ALIGNMENT_CHECK(alignment);
  // Prepare read
  gt_alignment_clear(alignment);
  // Parse FASTA/FASTQ record
  return gt_ifp_parse_fasta_fastq_read(buffered_fasta_input,
      alignment->tag,alignment->read,alignment->qualities,alignment->attributes);
}
GT_INLINE gt_status gt_input_fasta_parser_get_template(
    gt_buffered_input_file* const buffered_fasta_input,gt_template* const template,const bool paired_read) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  GT_TEMPLATE_CHECK(template);
  // Prepare read
  gt_template_clear(template,true);
  // Parse FASTA/FASTQ record (end/1)
  gt_status error_code;
  gt_alignment* alignment = gt_template_get_block_dyn(template,0);
  if ((error_code=gt_input_fasta_parser_get_alignment(buffered_fasta_input,alignment))!=GT_IFP_OK) {
    return error_code;
  }
  if (paired_read) {
    gt_alignment* alignment = gt_template_get_block_dyn(template,1);
    if (gt_buffered_input_file_eob(buffered_fasta_input)) return GT_IFP_FAIL;
    // Parse FASTA/FASTQ record (end/2)
    if ((error_code=gt_input_fasta_parser_get_alignment(buffered_fasta_input,alignment)) != GT_IFP_OK) {
      return error_code;
    }
    // Set pair attribute
    gt_attributes_add(template->attributes,GT_ATTR_ID_TAG_PAIR,
        gt_attributes_get(alignment->attributes,GT_ATTR_ID_TAG_PAIR),int64_t);
  } else {
    // Copy pairing information to tag attributes
    gt_attributes_copy(template->attributes,alignment->attributes); // FIXME: Faster...
    int64_t pair = GT_PAIR_BOTH;
    gt_attributes_add(template->attributes,GT_ATTR_ID_TAG_PAIR,&pair,int64_t);
  }
  return GT_IFP_OK;
}

#define GT_INPUT_MULTIFASTA_RETURN_ERROR(error_code) \
  gt_string_delete(buffer);  \
  gt_segmented_sequence_delete(seg_seq); \
  return GT_IFP_PE_TAG_BAD_BEGINNING

GT_INLINE gt_status gt_input_multifasta_parser_get_archive(
    gt_input_file* const input_multifasta_file,gt_sequence_archive* const sequence_archive) {
  GT_INPUT_FILE_CHECK(input_multifasta_file);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Clear sequence archive (Why? Could be done in an accumulative fashion)
  // gt_sequence_archive_clear(sequence_archive);
  // Check the file. Reload buffer if needed
  if (gt_input_file_eof(input_multifasta_file)) return GT_IFP_OK;
  GT_INPUT_FILE_CHECK(input_multifasta_file);
  // Read all sequences
  gt_string* const buffer = gt_string_new(200); // TODO: Should be done in terms of dna_string
  while (!gt_input_file_eof(input_multifasta_file)) {
    gt_segmented_sequence* seg_seq = gt_segmented_sequence_new();
    // Parse TAG
    if (gt_input_file_get_current_char(input_multifasta_file)!=GT_IFP_FASTA_TAG_BEGIN) {
      GT_INPUT_MULTIFASTA_RETURN_ERROR(GT_IFP_PE_TAG_BAD_BEGINNING);
    }
    gt_input_file_next_char(input_multifasta_file);
    while (gt_expect_true(!gt_input_file_eof(input_multifasta_file) &&
        gt_input_file_get_current_char(input_multifasta_file)!=EOL &&
        gt_input_file_get_current_char(input_multifasta_file)!=DOS_EOL)) {
      gt_string_append_char(seg_seq->seq_name,gt_input_file_get_current_char(input_multifasta_file));
      gt_input_file_next_char(input_multifasta_file);
    }
    gt_input_file_skip_eol(input_multifasta_file,NULL);
    gt_string_append_eos(seg_seq->seq_name);
    // Parse all the sequence lines
    while (gt_expect_true(!gt_input_file_eof(input_multifasta_file))) {
      if (gt_expect_true(gt_is_iupac_code(gt_input_file_get_current_char(input_multifasta_file)))) {
        // TODO: gt_get_dna_normalized should be transparent
        gt_string_append_char(buffer,
            gt_get_dna_normalized(gt_input_file_get_current_char(input_multifasta_file)));
        gt_input_file_next_char(input_multifasta_file);
      } else if (gt_input_file_get_current_char(input_multifasta_file)==EOL ||
                 gt_input_file_get_current_char(input_multifasta_file)==DOS_EOL) {
        // Append string and clear buffer (next line!)
        gt_segmented_sequence_append_string(seg_seq,
            gt_string_get_string(buffer),gt_string_get_length(buffer));
        gt_string_clear(buffer);
        gt_input_file_skip_eol(input_multifasta_file,NULL);
      } else if (gt_input_file_get_current_char(input_multifasta_file)==GT_IFP_FASTA_TAG_BEGIN) {
        break; // Next sequence !
      }
    }
    // Store the parsed sequence
    gt_sequence_archive_add_segmented_sequence(sequence_archive,seg_seq);
  }
  gt_string_delete(buffer);
  return GT_IFP_OK;
}
/*
 * Synch read of blocks
 */
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks(
    gt_buffered_input_file* const buffered_input1,gt_buffered_input_file* const buffered_input2,pthread_mutex_t* const input_mutex) {
  return gt_input_fasta_parser_synch_blocks_va(input_mutex,2,buffered_input1,buffered_input2);
}
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks_v(
    pthread_mutex_t* const input_mutex,uint64_t num_inputs,gt_buffered_input_file* const buffered_input,va_list v_args) {
  GT_NULL_CHECK(input_mutex);
  GT_ZERO_CHECK(num_inputs);
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input);
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed (synch)
  if (gt_buffered_input_file_eob(buffered_input)) {
    // Dump buffer if BOF it attached to the input, and get new out block (always FIRST)
    gt_buffered_input_file_dump_attached_buffers(buffered_input->attached_buffered_output_file);
    // Read synch blocks & Reload all the 'buffered_fasta_input' files
    GT_BEGIN_MUTEX_SECTION(*input_mutex) {
      if ((error_code=gt_input_fasta_parser_reload_buffer(buffered_input))!=GT_IFP_OK) {
        GT_END_MUTEX_SECTION(*input_mutex);
        return error_code;
      }
      while ((--num_inputs) > 0) {
        gt_buffered_input_file* remaining_buffered_input = va_arg(v_args,gt_buffered_input_file*);
        GT_BUFFERED_INPUT_FILE_CHECK(remaining_buffered_input);
        if ((error_code=gt_input_fasta_parser_reload_buffer(remaining_buffered_input))!=GT_IFP_OK) {
          GT_END_MUTEX_SECTION(*input_mutex);
          return error_code;
        }
      }
    } GT_END_MUTEX_SECTION(*input_mutex);
  }
  return GT_IFP_OK;
}
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,const uint64_t num_inputs,gt_buffered_input_file* const buffered_input,...) {
  GT_NULL_CHECK(input_mutex);
  GT_ZERO_CHECK(num_inputs);
  va_list v_args;
  va_start(v_args,buffered_input);
  gt_status error_code = gt_input_fasta_parser_synch_blocks_v(input_mutex,num_inputs,buffered_input,v_args);
  va_end(v_args);
  return error_code;
}
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks_a(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file** const buffered_input,const uint64_t num_inputs) {
  GT_NULL_CHECK(input_mutex);
  GT_NULL_CHECK(buffered_input);
  GT_ZERO_CHECK(num_inputs);
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed (synch)
  if (gt_buffered_input_file_eob(buffered_input[0])) {
    // Dump buffer if BOF it attached to the input, and get new out block (always FIRST)
    GT_BUFFERED_INPUT_FILE_CHECK(buffered_input[0]);
    gt_buffered_input_file_dump_attached_buffers(buffered_input[0]->attached_buffered_output_file);
    // Read synch blocks & Reload all the 'buffered_input' files
    GT_BEGIN_MUTEX_SECTION(*input_mutex) {
      uint64_t i;
      for (i=0;i<num_inputs;++i) {
        GT_BUFFERED_INPUT_FILE_CHECK(buffered_input[i]);
        if ((error_code=gt_input_fasta_parser_reload_buffer(buffered_input[i]))!=GT_IFP_OK) {
          GT_END_MUTEX_SECTION(*input_mutex);
          return error_code;
        }
      }
    } GT_END_MUTEX_SECTION(*input_mutex);
  }
  return GT_IFP_OK;
}
