/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_map_parser.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_input_map_parser.h"

#define GT_IMP_NUM_LINES          GT_NUM_LINES_10K
#define GT_IMP_SUBSET_NUM_LINES   GT_NUM_LINES_10K
#define GT_IMP_NUM_INIT_TAG_CHARS 200

/*
 * Useful macros
 */
// MAP/MAPQ related
#define gt_is_valid_template_separator(character) \
  (character==' ')
#define gt_is_valid_counter_separator(character) \
  (character=='+' || character==':' || character=='x')

/*
 * Map Parser Attributes
 */
GT_INLINE gt_map_parser_attributes* gt_input_map_parser_attributes_new(const bool force_read_paired) {
  gt_map_parser_attributes* attributes = gt_alloc(gt_map_parser_attributes);
  gt_input_map_parser_attributes_reset_defaults(attributes);
  gt_input_map_parser_attributes_set_paired(attributes,force_read_paired);
  return attributes;
}
GT_INLINE void gt_input_map_parser_attributes_delete(gt_map_parser_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  gt_free(attributes);
}
GT_INLINE void gt_input_map_parser_attributes_reset_defaults(gt_map_parser_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->max_parsed_maps = GT_ALL;
  attributes->force_read_paired = false;
  attributes->src_text = NULL;
  attributes->skip_based_model=false;
  attributes->remove_duplicates=false;
}
GT_INLINE bool gt_input_map_parser_attributes_is_paired(gt_map_parser_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  return attributes->force_read_paired;
}
GT_INLINE void gt_input_map_parser_attributes_set_paired(gt_map_parser_attributes* const attributes,const bool force_read_paired) {
  GT_NULL_CHECK(attributes);
  attributes->force_read_paired = force_read_paired;
}
GT_INLINE void gt_input_map_parser_attributes_set_max_parsed_maps(gt_map_parser_attributes* const attributes,const uint64_t max_parsed_maps) {
  GT_NULL_CHECK(attributes);
  attributes->max_parsed_maps = max_parsed_maps;
}
GT_INLINE void gt_input_map_parser_attributes_set_src_text(gt_map_parser_attributes* const attributes,gt_string* const src_text) {
  GT_NULL_CHECK(attributes);
  attributes->src_text = src_text;
}
GT_INLINE void gt_input_map_parser_attributes_set_skip_model(gt_map_parser_attributes* const attributes,const bool skip_based_model) {
  GT_NULL_CHECK(attributes);
  attributes->skip_based_model = skip_based_model;
}
GT_INLINE void gt_input_map_parser_attributes_set_duplicates_removal(gt_map_parser_attributes* const attributes,const bool remove_duplicates) {
  GT_NULL_CHECK(attributes);
  attributes->remove_duplicates = remove_duplicates;
}

/*
 * MAP File Format test
 */
GT_INLINE bool gt_input_map_parser_test_map(
    char* const file_name,const uint64_t line_num,char* const buffer,const uint64_t buffer_size,
    gt_map_file_format* const map_file_format,const bool show_errors) {
  // Count tabs
  uint64_t buffer_pos=0, num_tabs=0;
  int64_t begin_f2=-1, end_f2=-1;
  int64_t begin_f3=-1, end_f3=-1;
  int64_t begin_f4=-1, end_f4=-1;
  while (buffer_pos<buffer_size) {
    const char c = buffer[buffer_pos];
    // Check TAB/EOL
    if (gt_expect_false(c==TAB)) {
      if (num_tabs==0) {begin_f2=buffer_pos+1;}
      if (num_tabs==1) {end_f2=buffer_pos; begin_f3=buffer_pos+1;}
      if (num_tabs==2) {end_f3=buffer_pos; begin_f4=buffer_pos+1;}
      if (num_tabs==3) {end_f4=buffer_pos;}
      ++num_tabs;
    } else if (gt_expect_false(c==EOL)) {
      break;
    }
    ++buffer_pos;
  }
  // Set MAP type {MAP,MAPQ}
  const bool contains_qualities = (num_tabs==4);
  const bool null_qualities = (end_f3-begin_f3)==0;
  // Check MAP file
  //   MAP:   TAG\tREAD\tCOUNTERS\tMAPS
  //   MAPQ:  TAG\tREAD\tQUALS\tCOUNTERS\tMAPS
  // Required conditions:
  //   (1) 3|4 TABS
  //   (2) length(read)==length(quals)
  if (num_tabs!=3 && num_tabs!=4) {
    gt_cond_error(show_errors,PARSE_MAP_BAD_NUMBER_FIELDS,file_name,line_num,num_tabs);
    return false;
  } else if (gt_expect_false(num_tabs==4 && !null_qualities && (end_f2-begin_f2)!=(end_f3-begin_f3))) {
    gt_cond_error(show_errors,PARSE_MAP_BAD_READ_QUAL_LENGTH,
        file_name,line_num,end_f2-begin_f2,end_f3-begin_f3);
    return false;
  }
  // Check counters format
  const uint64_t begin_counters = (!contains_qualities) ? begin_f3 : begin_f4;
  const uint64_t end_counters = (!contains_qualities) ? end_f3 : end_f4;
  bool is_prev_number = false;
  for (buffer_pos=begin_counters;buffer_pos<end_counters;++buffer_pos) {
    const char c = buffer[buffer_pos];
    if (gt_is_number(c)) {
      is_prev_number = true;
    } else {
      if (gt_expect_false(!is_prev_number || !gt_is_valid_counter_separator(c)) ) {
        gt_cond_error(show_errors,PARSE_MAP_COUNTERS,file_name,line_num,buffer_pos+1);
        return false;
      }
      is_prev_number = false;
    }
  }
  // Extract extra MAP format information
  uint64_t num_blocks=1;
  char block_separator = 0;
  // Detect read template blocks
  const uint64_t f3_f4_diff = (begin_f4-begin_f3);
  for (buffer_pos=begin_f2;buffer_pos<end_f2;++buffer_pos) {
    const char c = buffer[buffer_pos];
    if (!gt_is_dna(c)) { // Check template block separator
      if (gt_expect_false(!gt_is_valid_template_separator(c))) {
        gt_cond_error(show_errors,PARSE_MAP_BAD_TEMPLATE_SEP,
            file_name,line_num,buffer_pos+1,c,"Not a valid separator");
        return false;
      } else if (gt_expect_false(block_separator!=0 && c!=block_separator)) {
        gt_cond_error(show_errors,PARSE_MAP_BAD_TEMPLATE_SEP,
            file_name,line_num,buffer_pos+1,c,"Not consistent with previous separator");
        return false;
      } else if (gt_expect_false(contains_qualities && !null_qualities && c!=buffer[buffer_pos+f3_f4_diff])) {
        gt_cond_error(show_errors,PARSE_MAP_BAD_TEMPLATE_SEP,
            file_name,line_num,buffer_pos+1,c,"Not synchronized with qualities separator");
        return false;
      }
      block_separator = c;
      ++num_blocks;
    }
  }
  // Set extra format attributes
  map_file_format->contains_qualities = contains_qualities;
  // TODO: IF qualities check (num_blocks(read==qualities)) ... and use DIFF_TEMPLATE_BLOCKS
  return true;
}
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors) {
  return gt_input_map_parser_test_map(
      gt_input_file_get_file_name(input_file),input_file->processed_lines+1,(char*)input_file->file_buffer,
      input_file->buffer_size,map_file_format,show_errors);
}
GT_INLINE gt_status gt_input_map_parser_check_map_file_format(gt_buffered_input_file* const buffered_input_file) {
  gt_input_file* const input_file = buffered_input_file->input_file;
  if (gt_expect_false(input_file->file_format==FILE_FORMAT_UNKNOWN)) { // Unknown
    gt_map_file_format map_type;
    if (!gt_input_map_parser_test_map(
        gt_input_file_get_file_name(input_file),buffered_input_file->current_line_num,
        gt_vector_get_mem(buffered_input_file->block_buffer,char),
        gt_vector_get_used(buffered_input_file->block_buffer),&map_type,true)) {
      return GT_IMP_PE_WRONG_FILE_FORMAT;
    }
    input_file->file_format = MAP;
    input_file->map_type = map_type;
  } else if (gt_expect_false(input_file->file_format!=MAP)) {
    return GT_IMP_PE_WRONG_FILE_FORMAT;
  }
  return 0;
}

/*
 * MAP File basics
 */
/* Error handler */
GT_INLINE void gt_input_map_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code) {
  // Display textual error msg
  const char* const file_name = (buffered_map_input != NULL) ?
      gt_input_file_get_file_name(buffered_map_input->input_file) : "<<LazyParsing>>";
  if ((buffered_map_input == NULL)) {
    line_num = 0; column_pos = 0;
  }
  switch (error_code) {
    case 0: /* No error */ break;
    case GT_IMP_PE_WRONG_FILE_FORMAT: gt_error(PARSE_MAP_BAD_FILE_FORMAT,file_name,line_num); break;
    case GT_IMP_PE_NOT_IMPLEMENTED: gt_error(PARSE_MAP_NOT_IMPLEMENTED,file_name,line_num,column_pos); break;
    case GT_IMP_PE_PREMATURE_EOL: gt_error(PARSE_MAP_PREMATURE_EOL,file_name,line_num,column_pos); break;
    case GT_IMP_PE_BAD_NUMBER_OF_BLOCKS: /* (blocks(read)!=blocks(quals) || parse_alignment=>(num_blocks==1)) */
      gt_error(PARSE_MAP_BAD_NUMBER_OF_BLOCKS,file_name,line_num,column_pos);
      break;
    case GT_IMP_PE_BAD_CHARACTER: gt_error(PARSE_MAP_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_READ_BAD_CHARACTER: gt_error(PARSE_MAP_READ_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_QUAL_BAD_SEPARATOR: gt_error(PARSE_MAP_QUAL_BAD_SEPARATOR,file_name,line_num,column_pos); break;
    case GT_IMP_PE_QUAL_BAD_LENGTH: gt_error(PARSE_MAP_QUAL_BAD_LENGTH,file_name,line_num,column_pos); break;
    case GT_IMP_PE_QUAL_BAD_CHARACTER: gt_error(PARSE_MAP_QUAL_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_COUNTERS_BAD_CHARACTER: gt_error(PARSE_MAP_COUNTERS_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS: gt_error(PARSE_MAP_BAD_NUMBER_OF_BLOCKS,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MAP_INCONSISTENT_BLOCKS: gt_error(PARSE_MAP_INCONSISTENT_BLOCKS,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MAP_BAD_CHARACTER: gt_error(PARSE_MAP_BAD_CHARACTER,file_name,line_num,column_pos); break;
//    case GT_IMP_PE_SPLIT_MAP_BAD_NUM_ACCEPTORS: gt_error(PARSE_MAP_SPLIT_MAP_BAD_NUM_ACCEPTORS,file_name,line_num,column_pos); break;
//    case GT_IMP_PE_SPLIT_MAP_BAD_NUM_DONORS: gt_error(PARSE_MAP_SPLIT_MAP_BAD_NUM_DONORS,file_name,line_num,column_pos); break;
//    case GT_IMP_PE_MISMS_ALREADY_PARSED: gt_error(PARSE_MAP_MISMS_ALREADY_PARSED,file_name,line_num); break;
    case GT_IMP_PE_MISMS_BAD_CHARACTER: gt_error(PARSE_MAP_MISMS_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MISMS_BAD_MISMS_POS: gt_error(PARSE_MAP_MISMS_BAD_MISMS_POS,file_name,line_num,column_pos); break;
    default:
      gt_error(PARSE_MAP,file_name,line_num);
      break;
  }
}
/* MAP file. Skip record */
GT_INLINE void gt_input_map_parser_next_record(gt_buffered_input_file* const buffered_map_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  gt_buffered_input_file_skip_line(buffered_map_input);
}
/* Read last record's tag (for block synchronization purposes ) */
GT_INLINE gt_status gt_imp_get_tag_last_read(gt_buffered_input_file* const buffered_map_input,gt_string* const last_tag) {
  int64_t position = gt_vector_get_used(buffered_map_input->block_buffer)-1;
  char* text_line = gt_vector_get_last_elm(buffered_map_input->block_buffer,char);
  if (*text_line!=EOL) return GT_INPUT_STATUS_FAIL;
  // Skip EOL
  while (*text_line==EOL) {--text_line; --position;}
  // Skip Text and get previos EOL and TABS
  uint64_t last_tab = 0;
  while (position>0 && *text_line!=EOL) {
    if (*text_line==TAB) last_tab = position;
    --text_line; --position;
  }
  // Check end cases
  if (last_tab==0 || last_tab<=position) return GT_INPUT_STATUS_FAIL; // Sth is wrong
  // Set last record's tag
  if (*text_line==EOL) { ++text_line; ++position; }
  gt_string_set_nstring(last_tag,text_line,last_tab-position);
  return GT_INPUT_STATUS_OK;
}
GT_INLINE gt_status gt_imp_reload_buffer_matching_tag(
    gt_buffered_input_file* const buffered_map_input,gt_string* const reference_tag,const bool read_paired) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  // Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
  gt_buffered_input_file_dump_attached_buffers(buffered_map_input->attached_buffered_output_file);
  // Read new input block
  gt_input_file* const input_file = buffered_map_input->input_file;
  // Read lines
  if (gt_input_file_eof(input_file)) return GT_INPUT_STATUS_EOF;
  gt_input_file_lock(input_file);
  if (gt_input_file_eof(input_file)) {
    gt_input_file_unlock(input_file);
    return GT_INPUT_STATUS_EOF;
  }
  buffered_map_input->block_id = gt_input_file_get_next_id(input_file);
  buffered_map_input->current_line_num = input_file->processed_lines+1;
  gt_vector_clear(buffered_map_input->block_buffer); // Clear dst buffer
  // Read lines
  if (read_paired) gt_input_parse_tag_chomp_pairend_info(reference_tag);
  gt_string* const last_tag = gt_string_new(0);
  uint64_t lines_read, total_lines_read = 0;
  uint64_t num_blocks = 0, num_tabs = 0;
  do {
    if ((lines_read=gt_input_file_next_record(input_file,
        buffered_map_input->block_buffer,last_tag,&num_blocks,&num_tabs))==0) break;
    if (read_paired) gt_input_parse_tag_chomp_pairend_info(last_tag);
    ++total_lines_read;
  } while (!gt_string_equals(reference_tag,last_tag));
  if (read_paired && lines_read>0 && num_blocks%2!=0) { // Check paired read
    gt_input_file_next_line(input_file,buffered_map_input->block_buffer);
    ++total_lines_read;
  }
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffered_map_input->block_buffer);
  if (total_lines_read > 0 && *gt_vector_get_last_elm(buffered_map_input->block_buffer,char) != EOL) {
    gt_vector_insert(buffered_map_input->block_buffer,EOL,char);
  }
  input_file->processed_lines+=total_lines_read;
  buffered_map_input->lines_in_buffer = total_lines_read;
  gt_input_file_unlock(input_file);
  // Setup the block
  buffered_map_input->cursor = gt_vector_get_mem(buffered_map_input->block_buffer,char);
  // Assign block ID
  gt_buffered_input_file_set_id_attached_buffers(buffered_map_input->attached_buffered_output_file,buffered_map_input->block_id);
  // Free
  gt_string_delete(last_tag);
  return GT_INPUT_STATUS_OK;
}
/* MAP file. Synchronized get block wrt to paired map records */
GT_INLINE gt_status gt_imp_get_block(
    gt_buffered_input_file* const buffered_map_input,const uint64_t num_records) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  gt_input_file* const input_file = buffered_map_input->input_file;
  // Read lines
  if (gt_input_file_eof(input_file)) return GT_INPUT_STATUS_EOF;
  gt_input_file_lock(input_file);
  if (gt_input_file_eof(input_file)) {
    gt_input_file_unlock(input_file);
    return GT_INPUT_STATUS_EOF;
  }
  buffered_map_input->block_id = gt_input_file_get_next_id(input_file);
  buffered_map_input->current_line_num = input_file->processed_lines+1;
  gt_vector_clear(buffered_map_input->block_buffer); // Clear dst buffer
  // Read lines
  uint64_t lines_read = 0, num_blocks = 0, num_tabs = 0;
  while ( (lines_read<num_records || num_blocks%2!=0) &&
      gt_input_file_next_record(input_file,buffered_map_input->block_buffer,NULL,&num_blocks,&num_tabs) ) ++lines_read;
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffered_map_input->block_buffer);
  if (lines_read > 0 && *gt_vector_get_last_elm(buffered_map_input->block_buffer,char) != EOL) {
    gt_vector_insert(buffered_map_input->block_buffer,EOL,char);
  }
  input_file->processed_lines+=lines_read;
  buffered_map_input->lines_in_buffer = lines_read;
  gt_input_file_unlock(input_file);
  // Setup the block
  buffered_map_input->cursor = gt_vector_get_mem(buffered_map_input->block_buffer,char);
  return buffered_map_input->lines_in_buffer;
}
/* MAP file. Reload internal buffer */
GT_INLINE gt_status gt_input_map_parser_reload_buffer(
    gt_buffered_input_file* const buffered_map_input,const bool synchronized_map,const uint64_t num_lines) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  // Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
  gt_buffered_input_file_dump_attached_buffers(buffered_map_input->attached_buffered_output_file);
  // Read new input block
  const uint64_t read_lines = (synchronized_map) ?
      gt_imp_get_block(buffered_map_input,num_lines):
      gt_buffered_input_file_get_block(buffered_map_input,num_lines);
  if (gt_expect_false(read_lines==0)) return GT_INPUT_STATUS_EOF;
  // Assign block ID
  gt_buffered_input_file_set_id_attached_buffers(buffered_map_input->attached_buffered_output_file,buffered_map_input->block_id);
  return GT_INPUT_STATUS_OK;
}

/*
 * MAP format. Basic building block for parsing
 */
GT_INLINE gt_status gt_imp_tag(const char** const text_line,gt_string* const tag,gt_attributes* const attributes) {
  // Delimit the tag
  gt_input_parse_tag(text_line,tag,attributes);
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  return 0;
}
GT_INLINE gt_status gt_imp_read_block(const char** const text_line,gt_string* const read_block) {
  // Read READ_BLOCK
  const char* const read_block_begin = *text_line;
  while (gt_expect_true(**text_line!=TAB && !gt_is_valid_template_separator(**text_line) && !GT_IS_EOL(text_line))) {
    if (gt_expect_false(!gt_is_dna(**text_line))) return GT_IMP_PE_READ_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
  }
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  // Copy string
  gt_string_set_nstring_static(read_block,read_block_begin,(*text_line-read_block_begin));
  // Place cursor at beginning of the next field
  gt_status return_status;
  if (**text_line==TAB) {
    return_status = GT_IMP_PE_EOB;
    GT_NEXT_CHAR(text_line);
  } else if (gt_is_valid_template_separator(**text_line)) {
    return_status = GT_IMP_PE_PENDING_BLOCKS;
    GT_NEXT_CHAR(text_line);
  } else {
    return_status = GT_IMP_PE_READ_BAD_CHARACTER;
  }
  return return_status;
}
GT_INLINE gt_status gt_imp_qualities_block(const char** const text_line,gt_string* const qualities_block) {
  // Read QUAL_BLOCK
  const char* const qualities_block_begin = *text_line;
  while (gt_expect_true(**text_line!=TAB && !gt_is_valid_template_separator(**text_line) && !GT_IS_EOL(text_line))) {
    if (gt_expect_false(!gt_is_valid_quality(**text_line))) return GT_IMP_PE_QUAL_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
  }
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  // Copy string
  gt_string_set_nstring_static(qualities_block,qualities_block_begin,(*text_line-qualities_block_begin));
  // Place cursor at beginning of the next field
  gt_status return_status;
  if (**text_line==TAB) {
    return_status = GT_IMP_PE_EOB;
    GT_NEXT_CHAR(text_line);
  } else if (gt_is_valid_template_separator(**text_line)) {
    return_status = GT_IMP_PE_PENDING_BLOCKS;
    GT_NEXT_CHAR(text_line);
  } else {
    return GT_IMP_PE_QUAL_BAD_SEPARATOR;
  }
  return return_status;
}
GT_INLINE gt_status gt_imp_counters(const char** const text_line,gt_vector* const counters,gt_attributes* const attributes) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_VECTOR_CHECK(counters);
  GT_ATTRIBUTES_CHECK(attributes);
  // Handle 'not-unique' flag
  if (**text_line==GT_MAP_COUNTS_NOT_UNIQUE) {
    bool not_unique = true;
    gt_attributes_add(attributes,GT_ATTR_ID_NOT_UNIQUE,&not_unique,bool);
    GT_NEXT_CHAR(text_line);
    if (gt_expect_false(**text_line!=TAB)) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
    return 0;
  }
  // Parse counters
  bool prev_char_was_sep = false, is_mcs_set = false;
  uint64_t number;
  while (gt_expect_true(**text_line!=TAB && !GT_IS_EOL(text_line))) {
    if (gt_is_number(**text_line)) {
      GT_PARSE_NUMBER(text_line,number);
      gt_vector_insert(counters,number,uint64_t);
      prev_char_was_sep = false;
    } else if (**text_line==GT_MAP_COUNTS_TIMES) { // 0x10:1:1
      if (prev_char_was_sep) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      if (gt_expect_false(gt_vector_get_used(counters)==0)) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      uint64_t multiplier, i;
      // Parse multiplier
      GT_NEXT_CHAR(text_line);
      if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
      if (!gt_is_number(**text_line)) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,multiplier);
      number = *gt_vector_get_last_elm(counters,uint64_t);
      for (i=1;i<multiplier;++i) {
        gt_vector_insert(counters,number,uint64_t);
      }
    } else if (**text_line==GT_MAP_MCS) {
      if (prev_char_was_sep || is_mcs_set) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      is_mcs_set = true;
      gt_attributes_add(attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA,&(gt_vector_get_used(counters)),uint64_t);
      GT_NEXT_CHAR(text_line);
      prev_char_was_sep = true;
    } else if (**text_line==GT_MAP_COUNTS_SEP) {
      if (prev_char_was_sep) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      if (*(*text_line+1)==GT_MAP_COUNTS_SEP) break;
      GT_NEXT_CHAR(text_line);
      prev_char_was_sep = true;
    } else {
      return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
    }
  }
  // Set default MCS
  if (!is_mcs_set) {
    gt_attributes_add(attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA,&(gt_vector_get_used(counters)),uint64_t);
  }
  // Parse attributes (if any)
  if (**text_line==GT_MAP_COUNTS_SEP && *(*text_line+1)==GT_MAP_COUNTS_SEP) { // 0:0::<Value>::<Value>
    return GT_IMP_PE_NOT_IMPLEMENTED; // Unlikely to be include into the grammar
  }
  return 0;
}
GT_INLINE gt_status gt_imp_parse_strand(const char** const text_line,gt_strand* const strand) {
  switch ((**text_line)) {
    case GT_MAP_STRAND_FORWARD_SYMBOL:
    case GT_MAP_STRAND_FORWARD_LETTER:
      *strand=FORWARD;
    break;
    case GT_MAP_STRAND_REVERSE_SYMBOL:
    case GT_MAP_STRAND_REVERSE_LETTER:
      *strand=REVERSE;
    break;
    default:
      return GT_IMP_PE_MAP_BAD_CHARACTER;
      break;
  }
  return 0;
}
GT_INLINE gt_status gt_imp_parse_mismatch_string_v0(const char** const text_line,gt_map* map,gt_map_parser_attributes* const map_parser_attr) {
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  /*
   * ReturnValues = { 0=OK, GT_IMP_PE_MISMS_BAD_CHARACTER, GT_IMP_PE_MISMS_BAD_MISMS_POS }
   * Parses GEM mismatch/CIGAR string v0
   *   OLD (v0): <+3>20A88C89C99
   */
  gt_map_clear_misms(map);
  // Store reference map (left-most position in the genome, strand, ...)
  gt_map* const head_map = map;
  const uint64_t global_length = gt_map_get_base_length(map);
  const uint64_t start_position = gt_map_get_position(map);
  const bool reverse_strand = (gt_map_get_strand(map)==REVERSE);
  // Auxiliary variables as to track the position in the read and the base length of the chunks
  uint64_t last_position = 0, last_cut_point = 0, num_blocks = 1;
  while ((**text_line)!=GT_MAP_NEXT && (**text_line)!=GT_MAP_SEP && !GT_IS_EOL(text_line) && (**text_line)!=GT_MAP_SCORE_GEMv0) {
    gt_misms misms;
    if (gt_is_dna((**text_line))) {
      /*
       * Mismatch
       */
      misms.misms_type = MISMS;
      misms.base = (**text_line);
      GT_NEXT_CHAR(text_line);
      // Parse Position
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,misms.position);
      if (gt_expect_false(misms.position<=last_position)) {
        return GT_IMP_PE_MISMS_BAD_MISMS_POS;
      }
      --misms.position; // Zero based position
      last_position = misms.position;
      misms.position -= last_cut_point; // Split-offset correction
      // Add Mismatch
      gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='<') {
      /*
       * PARSE Indel
       */
      bool is_splice;
      GT_NEXT_CHAR(text_line);
      // Parse operation [+-*]
      switch ((**text_line)) {
        case GT_MAP_INDEL_INSERTION:
          misms.misms_type = INS;
          is_splice = false;
          break;
        case GT_MAP_INDEL_DELETION:
          misms.misms_type = DEL;
          is_splice = false;
          break;
        case GT_MAP_INDEL_SPLICE:
          is_splice = true;
          break;
        default:
          return GT_IMP_PE_MISMS_BAD_CHARACTER;
          break;
      }
      GT_NEXT_CHAR(text_line);
      // Parse size
      uint64_t size, position;
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,size);
      // Parse Indel end ">"
      if (gt_expect_false((**text_line)!='>')) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_NEXT_CHAR(text_line);
      // Parse Position
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,position);
      if (gt_expect_false(position<=last_position)) {
        return GT_IMP_PE_MISMS_BAD_MISMS_POS;
      }
      --position; // Zero based position
      last_position = position;
      /*
       * Add Indel
       */
      if (gt_expect_true(!is_splice)) {
        misms.position = position-last_cut_point;
        misms.size = size;
        gt_map_add_misms(map,&misms);
      } else {
        // Close current map block
        gt_map_set_base_length(map,position-last_cut_point);
        last_cut_point = position;
        // Create a new map block
        gt_map* next_map = gt_map_new();
        gt_map_set_seq_name(next_map,gt_map_get_seq_name(map),gt_map_get_seq_name_length(map));
        gt_map_set_strand(next_map,gt_map_get_strand(map));
        gt_map_set_base_length(next_map,global_length-position);
        // Attach the next block
        gt_map_set_next_block(map,next_map,SPLICE,size); ++num_blocks;
        gt_map_set_position(next_map,gt_map_get_position(map)+gt_map_get_base_length(map)+size);
        // Swap maps & Reset length,position
        map = next_map;
      }
    } else { // ? Parsing error
      return GT_IMP_PE_MISMS_BAD_CHARACTER;
    }
  }
  // Adjust positions in case of SM in the reverse strand
  if (num_blocks>1 && reverse_strand) gt_map_reverse_blocks_positions(head_map,start_position);
  return 0;
}

GT_INLINE gt_status gt_imp_parse_mismatch_string_v1(const char** const text_line,gt_map* map,gt_map_parser_attributes* const map_parser_attr) {
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  /*
   * ReturnValues = { 0=OK, GT_IMP_PE_MISMS_BAD_CHARACTER }
   * Parses GEM mismatch/CIGAR string v1
   *   NEW(v1) ::= (5)43T46A9>24*
   *               33C9T30T24>1-(10)
   */
  gt_map_clear_misms(map);
  // Store reference map (left-most position in the genome, strand, ...)
  gt_map* const head_map = map;
  const uint64_t start_position = gt_map_get_position(map);
  const bool reverse_strand = (gt_map_get_strand(map)==REVERSE);
  // Auxiliary variables as to track the span of the map in the read and in the reference
  uint64_t read_span=0, reference_span=0, num_map_blocks=1;
  while ((**text_line)!=GT_MAP_NEXT && (**text_line)!=GT_MAP_SEP && !GT_IS_EOL(text_line)) {
    gt_misms misms;
    if (gt_is_number((**text_line))) {
      /*
       * Matching
       */
      uint64_t matching_characters;
      GT_PARSE_NUMBER(text_line,matching_characters);
      read_span+=matching_characters;
      reference_span+=matching_characters;
    } else if (gt_is_dna((**text_line))) {
      /*
       * Mismatch
       */
      misms.misms_type = MISMS;
      misms.base = (**text_line);
      misms.position = read_span;
      ++read_span; ++reference_span;
      GT_NEXT_CHAR(text_line);
      gt_map_add_misms(map,&misms); // Add Mismatch
    } else if ((**text_line)=='(') {
      /*
       * Trim
       */
      misms.misms_type = DEL;
      misms.position = read_span;
      GT_NEXT_CHAR(text_line);
      // Parse size
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,misms.size);
      read_span+=misms.size;
      // Parse Trim end ')'
      if (gt_expect_false((**text_line)!=')')) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_NEXT_CHAR(text_line);
      // Add Trim
      if (misms.size>0) gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='>') {
      /*
       * Indel/Skip
       */
      GT_NEXT_CHAR(text_line);
      // Parse size
      int64_t size;
      GT_PARSE_SIGNED_NUMBER_BLOCK(text_line,size) {
        if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
        GT_PARSE_NUMBER(text_line,size)
      } GT_PARSE_SIGNED_NUMBER_END_BLOCK(size);
      // Parse skip type
      if (!map_parser_attr->skip_based_model && (size > 0) &&
          ((**text_line)==GT_MAP_SKIP_POSITIVE || (**text_line)==GT_MAP_SKIP_NEGATIVE)) {
        /*
         * INS/DEL
         */
        misms.position = read_span;
        misms.size = size;
        if ((**text_line)==GT_MAP_SKIP_POSITIVE) {
          misms.misms_type = INS;
          reference_span+=misms.size;
        } else {
          misms.misms_type = DEL;
          read_span+=misms.size;
        }
        GT_NEXT_CHAR(text_line);
        // Add Indel/Skip
        gt_map_add_misms(map,&misms);
      } else {
        /*
         * NSKIP/SPLICE
         */
        gt_junction_t junction;
        switch ((**text_line)) {
          case GT_MAP_SKIP_POSITIVE: junction=POSITIVE_SKIP; break;
          case GT_MAP_SKIP_NEGATIVE: junction=NEGATIVE_SKIP; break;
          case GT_MAP_SKIP_SPLICE: junction=SPLICE; break;
          default: return GT_IMP_PE_MISMS_BAD_CHARACTER; break;
        }
        GT_NEXT_CHAR(text_line);
        // Create a new map block
        gt_map* const next_map = gt_map_new();
        gt_map_set_seq_name(next_map,gt_map_get_seq_name(map),gt_map_get_seq_name_length(map));
        gt_map_set_strand(next_map,gt_map_get_strand(map));
        // FIXME: gt_map_set_base_length(next_map,gt_map_get_base_length(map)-read_span);
        // Attach the next block & close current map block
        gt_map_set_base_length(map,read_span);
        gt_map_set_position(next_map,gt_map_get_position(map)+reference_span+size);
        gt_map_set_next_block(map,next_map,junction,size); ++num_map_blocks;
        // Swap maps & Reset reference_span,read_span
        map = next_map;
        read_span=0; reference_span=0;
      }
    } else { // ? Parsing error
     return GT_IMP_PE_MISMS_BAD_CHARACTER;
    }
  }
  // Set base length of the last map block
  gt_map_set_base_length(map,read_span);
  // Adjust positions in case of SM in the reverse strand
  if (num_map_blocks>1 && reverse_strand) gt_map_reverse_blocks_positions(head_map,start_position);
  return 0;
}
GT_INLINE gt_status gt_imp_parse_map_score_v0(const char** const text_line,uint64_t* const gt_score) {
  // Parse score1
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_PARSE_NUMBER(text_line,*gt_score);
  // Skip score2 (no use)
  if (gt_expect_false(**text_line!=GT_MAP_SCORE_SEP)) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_READ_UNTIL(text_line,!gt_is_number((**text_line)));
  return 0;
}
GT_INLINE gt_status gt_imp_parse_map_score_v1(const char** const text_line,uint64_t* const gt_score) {
  // Parse score (either a hex or decimal number)
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_PARSE_HEX_OR_DEC(text_line,*gt_score);
  return 0;
}
GT_INLINE gt_status gt_imp_next_element(const char** const text_line) {
  /*
   * ReturnValues = {
   *                 GT_IMP_PE_MAP_PENDING_MAPS, GT_IMP_PE_EOB, GT_IMP_PE_PENDING_BLOCKS
   *                 GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv0, GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv1,
   *                 GT_IMP_PE_MMAP_ATTRIBUTE_SCORE
   *                 GT_IMP_PE_MAP_BAD_CHARACTER
   *                 }
   */
  switch (**text_line) {
    case GT_MAP_NEXT: // ','
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MAP_PENDING_MAPS;
    case EOL:
    case EOS:
      return GT_IMP_PE_EOB; // '\n'
    case GT_MAP_SCORE_GEMv0: // '@10/1'
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv0;
    case GT_MAP_SEP:
      if ( (*(*text_line+1)) == GT_MAP_SEP) { // '::'
        if ( (*(*text_line+2)) == GT_MAP_SEP) { // ':::' (Attributes of the block group)
          if ( (*(*text_line+3)) != GT_MAP_SEP) {
            (*text_line)+=3;
            return GT_IMP_PE_MMAP_ATTRIBUTE_SCORE;
          } else {
            (*text_line)+=2;
            return GT_IMP_PE_PENDING_BLOCKS;
          }
        } else { // '::?'
          (*text_line)+=2;
          return GT_IMP_PE_PENDING_BLOCKS;
        }
      } else if (gt_is_number( (*(*text_line+1)) )) { // ':100'
        GT_NEXT_CHAR(text_line);
        return GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv1;
      } else { // ':?'
        return GT_IMP_PE_MAP_BAD_CHARACTER;
      }
    default:
      return GT_IMP_PE_MAP_BAD_CHARACTER;
  }
}
#define GT_IMP_PARSE_SPLIT_MAP_CLEAN1__RETURN(error_code) { gt_map_delete(donor_map); return error_code; }
#define GT_IMP_PARSE_SPLIT_MAP_CLEAN2__RETURN(error_code) { gt_map_delete(donor_map); gt_map_delete(acceptor_map); return error_code; }
#define GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line) ((**text_line)==GT_MAP_SPLITMAP_NEXT_GEMv0_0 || (**text_line)==GT_MAP_SPLITMAP_NEXT_GEMv0_1)
GT_INLINE gt_status gt_imp_parse_split_map_v0(const char** const text_line,gt_map** const split_map,const uint64_t read_base_length) {
  /*
   * ReturnValues = { GT_IMP_PE_MAP_BAD_CHARACTER, GT_IMP_PE_PREMATURE_EOL, OK=0 }
   */
  GT_NULL_CHECK(text_line);
  GT_NULL_CHECK(split_map);
  //
  // Split-map format (v0)
  //   [23]=chr6:R31322884~chr6:R31237276
  //   [70-71]=chr1:F188862944~chr19:F[53208292-53208293]
  //   [23-50]=chr1:F[188862944-188868041]~chr19:F53208292
  // But also we have ....
  //   [31;35]=chr16:R[2503415;2503411]~chr16:R2503271
  //   [30;34]=chr10:F74776624~chr10:F[74790025;74790029]
  // And also ...
  //   [25-26]=chr12:F95317898~chr12:F[95317924-95317925],[25-26]=chr11:R[100658012-100658011]~chr12:F[95317924-95317925]
  // Life is like a box of chocolates ...
  /*
   * Parse split-points
   */
  if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_OPEN_GEMv0)) return GT_IMP_PE_MAP_BAD_CHARACTER;
  // Create the SM
  gt_map* const donor_map = gt_map_new();
  // Read split-points
  uint64_t sm_position;
  bool sm_elm_parsed = false;
  do {
    GT_NEXT_CHAR(text_line);
    if (!gt_is_number((**text_line))) return GT_IMP_PE_MAP_BAD_CHARACTER;
    GT_PARSE_NUMBER(text_line,sm_position);
    if (!sm_elm_parsed) {
      gt_map_set_base_length(donor_map,sm_position);
      sm_elm_parsed=true;
    }
  } while (GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line));
  // Read closing split-points and definition symbol
  if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_CLOSE_GEMv0)) GT_IMP_PARSE_SPLIT_MAP_CLEAN1__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
  GT_NEXT_CHAR(text_line);
  if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_DEF_GEMv0)) GT_IMP_PARSE_SPLIT_MAP_CLEAN1__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
  GT_NEXT_CHAR(text_line);
  /*
   * Parse donor
   */
  // Read TAG
  const char* const donor_name = *text_line;
  GT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP);
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  gt_map_set_seq_name(donor_map,donor_name,(*text_line-donor_name));
  GT_NEXT_CHAR(text_line);
  // Read Strand
  gt_status error_code;
  if ((error_code=gt_imp_parse_strand(text_line,&(donor_map->strand)))) GT_IMP_PARSE_SPLIT_MAP_CLEAN1__RETURN(error_code);
  GT_NEXT_CHAR(text_line);
  // Detect multiple donor position
  bool multiple_sm;
  if ((**text_line)==GT_MAP_SPLITMAP_OPEN_GEMv0) {
    GT_NEXT_CHAR(text_line); // [31;35]=chr16:R[2503415;2503411]~chr16:R2503271
    multiple_sm = true;
  } else {
    multiple_sm = false;
  }
  // Read Position
  sm_elm_parsed = false;
  do {
    if (sm_elm_parsed) GT_NEXT_CHAR(text_line);
    if (gt_expect_false(!gt_is_number((**text_line)))) GT_IMP_PARSE_SPLIT_MAP_CLEAN1__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
    GT_PARSE_NUMBER(text_line,sm_position);
    if (!sm_elm_parsed) { // Set donor
      gt_map_set_position(donor_map,sm_position);
      sm_elm_parsed=true;
    }
  } while (GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line));
  // Close multiple donor syntax
  if (multiple_sm) {
    if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_CLOSE_GEMv0)) GT_IMP_PARSE_SPLIT_MAP_CLEAN1__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
    GT_NEXT_CHAR(text_line);
  }
  // Read separator (~)
  if (gt_expect_false(**text_line!=GT_MAP_SPLITMAP_SEP_GEMv0)) GT_IMP_PARSE_SPLIT_MAP_CLEAN1__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
  GT_NEXT_CHAR(text_line);
  /*
   * Parse acceptor(s)
   */
  // Read acceptor's TAG
  gt_map* const acceptor_map = gt_map_new();
  const char* const acceptor_name = *text_line;
  GT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP);
  if (GT_IS_EOL(text_line)) GT_IMP_PARSE_SPLIT_MAP_CLEAN2__RETURN(GT_IMP_PE_PREMATURE_EOL);
  gt_map_set_seq_name(acceptor_map,acceptor_name,(*text_line-acceptor_name));
  GT_NEXT_CHAR(text_line);
  // Read acceptor's Strand
  if ((error_code=gt_imp_parse_strand(text_line,&(acceptor_map->strand)))) GT_IMP_PARSE_SPLIT_MAP_CLEAN2__RETURN(error_code);
  GT_NEXT_CHAR(text_line);
  // Detect multiple donor position
  if (gt_expect_false((**text_line)==GT_MAP_SPLITMAP_OPEN_GEMv0)) { // [30;34]=chr10:F74776624~chr10:F[74790025;74790029]
    GT_NEXT_CHAR(text_line);
    multiple_sm = true;
  } else {
    multiple_sm = false;
  }
  // Parse acceptor position
  sm_elm_parsed = false;
  do {
    // Parse acceptor position
    if (sm_elm_parsed) GT_NEXT_CHAR(text_line);
    if (gt_expect_false(!gt_is_number((**text_line)))) GT_IMP_PARSE_SPLIT_MAP_CLEAN2__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
    GT_PARSE_NUMBER(text_line,sm_position);
    if (!sm_elm_parsed) { // Store split-map's position
      gt_map_set_position(acceptor_map,sm_position);
      sm_elm_parsed=true;
    }
  } while (GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line));
  // Close multiple acceptor syntax
  if (multiple_sm) {
    if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_CLOSE_GEMv0)) GT_IMP_PARSE_SPLIT_MAP_CLEAN2__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
    GT_NEXT_CHAR(text_line);
  }
  // Link both donor & acceptor
  gt_map_set_base_length(acceptor_map,read_base_length-gt_map_get_base_length(donor_map));
//  if (gt_map_get_strand(donor_map)==FORWARD) {
    gt_map_set_next_block(donor_map,acceptor_map,SPLICE,(int64_t)gt_map_get_position(acceptor_map)-(int64_t)gt_map_get_position(donor_map)+(int64_t)gt_map_get_length(donor_map));
//  } else {
//    const uint64_t acceptor_map_length = (read_base_length!=UINT64_MAX) ? gt_map_get_length(acceptor_map) : 0;
//    gt_map_set_next_block(acceptor_map,donor_map,SPLICE,(int64_t)gt_map_get_position(donor_map)-(int64_t)gt_map_get_position(acceptor_map)+(int64_t)acceptor_map_length);
//  }
  // Add the SM
//  if (gt_map_get_strand(donor_map)==FORWARD) {
    *split_map = donor_map;
//  } else {
//    *split_map = acceptor_map;
//  }
  // Return
  return 0;
}
#define GT_IMP_PARSE_MAP_CLEAN__RETURN(error_code) { \
  gt_map_delete(map); return error_code; \
}
GT_INLINE gt_status gt_imp_parse_map(
    const char** const text_line,gt_map** const return_map,
    const uint64_t read_base_length,gt_map_parser_attributes* const map_parser_attr) {
  /*
   * ReturnValues = {
   *                 GT_IMP_PE_PENDING_BLOCKS, GT_IMP_PE_MAP_PENDING_MAPS, GT_IMP_PE_EOB,
   *                 GT_IMP_PE_MMAP_ATTRIBUTE_SCORE
   *                 GT_IMP_PE_PREMATURE_EOL, GT_IMP_PE_MAP_BAD_CHARACTER, GT_IMP_PE_MISMS_BAD_CHARACTER, GT_IMP_PE_MISMS_BAD_MISMS_POS
   *                }
   * Parse MAP:
   *   OLD(v0) ::= chr7:F127708134G27T88
   *   SplitMap_OLD(v0) ::= [70-71]=chr1:F188862944~chr19:F[53208292-53208293]
   *                        [30;34]=chr10:F74776624~chr10:F[74790025;74790029]
   *                        [25-26]=chr11:R[100658012-100658011]~chr12:F[95317924-95317925]
   *                        ...
   *   NEW(v2)={chr11:-:51590050:(5)43TTC5>5-46A9>24*}
   */
  GT_NULL_CHECK(text_line);
  GT_NULL_CHECK(map_parser_attr);
  // Check null map (Orphan/unpaired/errors/oddities/...)
  gt_status error_code;
  if (**text_line==GT_MAP_NEXT) { // Pending (Null map)
    GT_NEXT_CHAR(text_line);
    *return_map=NULL;
    return GT_IMP_PE_MAP_PENDING_MAPS;
  } else if (GT_IS_EOL(text_line)) { // Pending (Null map)
    *return_map=NULL;
    return GT_IMP_PE_EOB;
  } else if (gt_expect_false(**text_line==GT_MAP_SEP && *(*text_line+1)==GT_MAP_SEP)) { // Orphan?
    GT_NEXT_CHAR(text_line);
    GT_NEXT_CHAR(text_line);
    *return_map=NULL;
    if (**text_line!=GT_MAP_SEP) {
      return GT_IMP_PE_PENDING_BLOCKS;
    } else {
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MMAP_ATTRIBUTE_SCORE;
    }
  } else if (gt_expect_false((**text_line)==GT_MAP_SPLITMAP_OPEN_GEMv0)) { // Parse Old Split-Maps
    if ((error_code=gt_imp_parse_split_map_v0(text_line,return_map,read_base_length))) return error_code;
  } else {
    /*
     * Parse MAP (Regular Map... for whatever that means)
     */
    gt_map* const map = gt_map_new();
    gt_map_set_base_length(map,read_base_length); // Tentative base length (for GEMv0)
    // Read TAG
    const char* const seq_name_start = *text_line;
    GT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP);
    if (GT_IS_EOL(text_line)) GT_IMP_PARSE_MAP_CLEAN__RETURN(GT_IMP_PE_PREMATURE_EOL);
    gt_map_set_seq_name(map,seq_name_start,(*text_line-seq_name_start));
    GT_NEXT_CHAR(text_line); // Separator ':'
    // Read Strand
    if ((error_code=gt_imp_parse_strand(text_line,&(map->strand)))) GT_IMP_PARSE_MAP_CLEAN__RETURN(error_code);
    GT_NEXT_CHAR(text_line); // Separator ':'
    // Determine format version
    gt_misms_string_t misms_format;
    if ((**text_line)==GT_MAP_SEP) { // GEMv1
      misms_format = MISMATCH_STRING_GEMv1;
      GT_NEXT_CHAR(text_line); // Separator ':'
    } else if (gt_is_number((**text_line))) { // GEMv0
      misms_format = MISMATCH_STRING_GEMv0;
    } else { // ?
      return GT_IMP_PE_MAP_BAD_CHARACTER;
    }
    // Position
    if (gt_expect_false(!gt_is_number((**text_line)))) GT_IMP_PARSE_MAP_CLEAN__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
    GT_PARSE_NUMBER(text_line,map->position);
    // Synch with mismatch string (GEMv1)
    if (misms_format==MISMATCH_STRING_GEMv1) {
      if (gt_expect_false((**text_line)!=GT_MAP_SEP)) GT_IMP_PARSE_MAP_CLEAN__RETURN(GT_IMP_PE_MAP_BAD_CHARACTER);
      GT_NEXT_CHAR(text_line);
    }
    // Parse Mismatch String
    if (misms_format==MISMATCH_STRING_GEMv1) {
      error_code=gt_imp_parse_mismatch_string_v1(text_line,map,map_parser_attr);
    } else {
      error_code=gt_imp_parse_mismatch_string_v0(text_line,map,map_parser_attr);
    }
    if (error_code) GT_IMP_PARSE_MAP_CLEAN__RETURN(error_code);
    // Return the parsed maps
    *return_map=map;
  }
  /*
   * Parse map's attributes (if any) and spot next item.
   */
  bool parsed_score = false;
  while (true) {
    error_code = gt_imp_next_element(text_line);
    switch (error_code) {
      case GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv0:
        if (parsed_score) return GT_IMP_PE_MAP_BAD_CHARACTER;
        parsed_score = true;
        if ((error_code=gt_imp_parse_map_score_v0(text_line,&((*return_map)->gt_score) ))) return error_code;
        (*return_map)->phred_score = (*return_map)->gt_score;
        break;
      case GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv1:
        if (parsed_score) return GT_IMP_PE_MAP_BAD_CHARACTER;
        parsed_score = true;
        if ((error_code=gt_imp_parse_map_score_v1(text_line,&((*return_map)->gt_score) ))) return error_code;
        (*return_map)->phred_score = (*return_map)->gt_score;
        break;
      default:
        return error_code;
        break;
    }
  }
}
#define GT_IMP_PARSE_TEMPLATE_MAP_RESET(map_head,last_map_block,total_base_length) \
  map_head = NULL; last_map_block = NULL; total_base_length = 0;
#define GT_IMP_PARSE_TEMPLATE_MAP_NEXT_END(current_end_position) { \
  mmap[current_end_position] = map_head; /* Keep current end */ \
  ++current_end_position;  /* Inc end position */ \
  if (error_code==GT_IMP_PE_PENDING_BLOCKS && current_end_position>1) return GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS; \
  GT_IMP_PARSE_TEMPLATE_MAP_RESET(map_head,last_map_block,total_base_length); \
}
#define GT_IMP_PARSE_MAP_BEGIN_ERROR_HANDLER(error_code) \
  switch (error_code) { \
    case GT_IMP_PE_PENDING_BLOCKS: \
    case GT_IMP_PE_MAP_PENDING_MAPS: \
    case GT_IMP_PE_EOB: \
    case GT_IMP_PE_MMAP_ATTRIBUTE_SCORE: \
      break; \
    default:
#define GT_IMP_PARSE_MAP_END_ERROR_HANDLER(error_code) \
      return error_code; \
      break; \
  }
GT_INLINE gt_status gt_imp_parse_template_maps(
    const char** const text_line,gt_template* const template,gt_map_parser_attributes* const map_parser_attr) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK((*text_line));
  GT_TEMPLATE_CHECK(template);
  // Check null maps
  if ((**text_line)==GT_MAP_NONE) {
    GT_SKIP_LINE(text_line);
    return 0;
  }
  // Check max_num_maps to parse (stratum-wise)
  uint64_t max_num_maps = map_parser_attr->max_parsed_maps;
  if (max_num_maps<GT_ALL) {
    uint64_t strata;
    gt_counters_calculate_num_maps(gt_template_get_counters_vector(template),
        0,map_parser_attr->max_parsed_maps,&strata,&max_num_maps);
  }
  // Parse MAPS. Formats allowed:
  //   OLD (v0): chr7:F127708134G27T88::chr7:R127708509<+3>20A88C89C99
  //   NEW (v1): chr11:-:51590050:(5)43T46A9>24*::chr11:-:51579823:33C9T30T24>1-(10)
  //       (v2): chrX:-:53678:100::,::chr2:-:3445:100
  gt_status error_code = GT_IMP_PE_MAP_PENDING_MAPS;
  uint64_t num_maps_parsed = 0;
  uint64_t read_length[2] = {gt_string_get_length(gt_template_get_end1(template)->read),gt_string_get_length(gt_template_get_end2(template)->read)};
  while (error_code==GT_IMP_PE_MAP_PENDING_MAPS && num_maps_parsed<max_num_maps) {
    /*
     * Parse MMAP. Read all possible not-null blocks (quimeras/oddities/...)
     */
    gt_map *mmap[2];
    gt_map *map_head = NULL, *last_map_block = NULL, *map_parsed;
    uint64_t total_base_length = 0, current_end_position = 0;
    do {
      // Parse MapBlock
      map_parsed = NULL;
      error_code = gt_imp_parse_map(text_line,&map_parsed,read_length[current_end_position],map_parser_attr);
      // Check error_code
      switch (error_code) {
        case GT_IMP_PE_PENDING_BLOCKS:
        case GT_IMP_PE_MAP_PENDING_MAPS:
        case GT_IMP_PE_EOB:
        case GT_IMP_PE_MMAP_ATTRIBUTE_SCORE: break;
        default: // Sth went wrong
          if (map_parsed!=NULL) gt_map_delete(map_parsed);
          if (map_head!=NULL) gt_map_delete(map_head);
          return error_code;
      }
      if (map_parsed==NULL) {
        // Check blocks & base lengths
        if (total_base_length>0 && total_base_length!=read_length[current_end_position]) {
          return GT_IMP_PE_MAP_INCONSISTENT_BLOCKS;
        }
        GT_IMP_PARSE_TEMPLATE_MAP_NEXT_END(current_end_position);
      } else {
        if (map_head==NULL) {
          map_head = map_parsed;
          last_map_block = map_head;
        } else {
          gt_map* const last_added = map_parsed;
          gt_map_set_next_block(last_map_block,last_added,QUIMERA,0);
          last_map_block = last_added;
        }
        // Check blocks & base lengths
        total_base_length += gt_map_get_global_base_length(map_parsed);
        if (total_base_length > read_length[current_end_position]) {
          return GT_IMP_PE_MAP_INCONSISTENT_BLOCKS;
        }
        else if (total_base_length == read_length[current_end_position]) {
          GT_IMP_PARSE_TEMPLATE_MAP_NEXT_END(current_end_position);
        }
      }
    } while (error_code==GT_IMP_PE_PENDING_BLOCKS);
    if (current_end_position<2 && total_base_length != read_length[current_end_position]) {
      return GT_IMP_PE_MAP_INCONSISTENT_BLOCKS;
    }
    if (current_end_position!=2) {
      return GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS;
    }
    /*
     * Parse attributes (if any) and calculate template attributes
     */
    gt_mmap_attributes mmap_attr;
    gt_template_mmap_attributes_clear(&mmap_attr);
    mmap_attr.distance = gt_mmap_get_global_distance(mmap,2);
    if (error_code==GT_IMP_PE_MMAP_ATTRIBUTE_SCORE) {
      if ((error_code=gt_imp_parse_map_score_v1(text_line,&mmap_attr.gt_score))) return error_code;
      mmap_attr.phred_score = mmap_attr.gt_score;
      error_code = gt_imp_next_element(text_line);
    } else {
      mmap_attr.gt_score = GT_MAP_NO_GT_SCORE;
      mmap_attr.phred_score = 0;
    }
    /*
     * Add the MMap to the template and Maps to individual alignments
     *   FIXME: Should aim to erase duplicates ... Not only to strictly parse it
     */
    if (mmap[0]!=NULL) gt_alignment_add_map(gt_template_get_end1(template),mmap[0]);
    if (mmap[1]!=NULL) gt_alignment_add_map(gt_template_get_end2(template),mmap[1]);
    gt_template_add_mmap_array(template,mmap,&mmap_attr);
    ++num_maps_parsed;
    /*
     * Check error code
     */
    switch (error_code) {
      case GT_IMP_PE_EOB:
      case GT_IMP_PE_MAP_PENDING_MAPS:
        break;
      case GT_IMP_PE_PENDING_BLOCKS:
        if (mmap[0]!=NULL) gt_map_delete(mmap[0]);
        if (mmap[1]!=NULL) gt_map_delete(mmap[1]);
        return GT_IMP_PE_MAP_BAD_CHARACTER;
      default: // Sth went wrong
        if (mmap[0]!=NULL) gt_map_delete(mmap[0]);
        if (mmap[1]!=NULL) gt_map_delete(mmap[1]);
        return error_code;
    }
  }
  return 0;
}
GT_INLINE gt_status gt_imp_parse_alignment_maps(const char** const text_line,gt_alignment* alignment,gt_map_parser_attributes* const map_parser_attr) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK((*text_line));
  GT_ALIGNMENT_CHECK(alignment);
  // Check null maps
  if ((**text_line)==GT_MAP_NONE) {
    GT_SKIP_LINE(text_line);
    return 0;
  }
  // Check max_num_maps to parse (stratum-wise)
  uint64_t max_num_maps = map_parser_attr->max_parsed_maps;
  if (max_num_maps<GT_ALL) {
    uint64_t strata;
    gt_counters_calculate_num_maps(gt_alignment_get_counters_vector(alignment),
        0,map_parser_attr->max_parsed_maps,&strata,&max_num_maps);
  }
  // Parse MAPS. Formats allowed:
  //   OLD (v0): chr7:F127708134G27T88
  //   NEW (v1): chr11:-:51590050:(5)43T46A9>24*
  const uint64_t alignment_base_length = gt_alignment_get_read_length(alignment);
  uint64_t num_maps_parsed = 0;
  gt_status error_code = GT_IMP_PE_MAP_PENDING_MAPS;
  while (error_code==GT_IMP_PE_MAP_PENDING_MAPS && num_maps_parsed<max_num_maps) {
    /*
     * Parse MAP. Read all possible not-null blocks (quimeras/oddities/...)
     */
    gt_map *map_head = NULL, *last_map_block = NULL, *map_parsed;
    uint64_t total_base_length = 0;
    do {
      // Parse MapBlock
      map_parsed = NULL;
      error_code = gt_imp_parse_map(text_line,&map_parsed,alignment_base_length,map_parser_attr);
      if (map_parsed==NULL) return GT_IMP_PE_MAP_BAD_CHARACTER;
      // Check error_code
      switch (error_code) {
        case GT_IMP_PE_PENDING_BLOCKS:
        case GT_IMP_PE_MAP_PENDING_MAPS:
        case GT_IMP_PE_EOB:
        case GT_IMP_PE_MMAP_ATTRIBUTE_SCORE: break;
        default: // Sth went wrong
          if (map_parsed!=NULL) gt_map_delete(map_parsed);
          if (map_head!=NULL) gt_map_delete(map_head);
          return error_code;
      }
      // Check blocks & base lengths
      if (map_head==NULL) {
        map_head = map_parsed;
        last_map_block = map_parsed;
      } else {
        gt_map* const last_added = map_parsed;
        gt_map_set_next_block(last_map_block,last_added,QUIMERA,0);
        last_map_block = last_added;
      }
      total_base_length += gt_map_get_global_base_length(last_map_block);
      if (total_base_length > alignment_base_length) {
        return GT_IMP_PE_MAP_INCONSISTENT_BLOCKS;
      }
    } while (error_code==GT_IMP_PE_PENDING_BLOCKS);
    if (total_base_length != alignment_base_length) {
      return GT_IMP_PE_MAP_INCONSISTENT_BLOCKS;
    }
    /*
     * Store map into alignment. FIXME: Should aim to erase duplicates ... Not only to strictly parse it
     */
    gt_alignment_add_map(alignment,map_head);
    ++num_maps_parsed;
    /*
     * Check error code
     */
    bool parse_attr = true, parsed_score = false;
    while (parse_attr) {
      switch (error_code) {
        case GT_IMP_PE_EOB:
        case GT_IMP_PE_MAP_PENDING_MAPS:
          parse_attr = false;
          break;
        case GT_IMP_PE_MMAP_ATTRIBUTE_SCORE:
          if (parsed_score) return GT_IMP_PE_MAP_BAD_CHARACTER;
          if ((error_code=gt_imp_parse_map_score_v1(text_line,&map_head->gt_score))) return error_code;
          map_head->phred_score = map_head->gt_score;
          parsed_score = true;
          error_code = gt_imp_next_element(text_line);
          break;
        case GT_IMP_PE_PENDING_BLOCKS: // No more pending blocks
          if (map_head!=NULL) gt_map_delete(map_head);
          return GT_IMP_PE_MAP_BAD_CHARACTER;
        default: // Sth went wrong
          if (map_head!=NULL) gt_map_delete(map_head);
          return error_code;
      }
    }
  }
  return 0;
}
GT_INLINE gt_status gt_imp_map_blocks(const char** const text_line,gt_map** const map,gt_map_parser_attributes* const map_parser_attr) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK((*text_line));
  GT_NULL_CHECK(map);
  // Check null maps
  if ((**text_line)==GT_MAP_NONE) {
    GT_NEXT_CHAR(text_line);
    return 0;
  }
  /*
   * Parse MAP. Read all possible not-null blocks (quimeras/oddities/...)
   */
  gt_status error_code;
  gt_map *map_head = NULL, *last_map_block = NULL, *map_parsed;
  do {
    // Parse MapBlock
    map_parsed = NULL;
    error_code = gt_imp_parse_map(text_line,&map_parsed,UINT32_MAX,map_parser_attr);
    if (map_parsed==NULL) return GT_IMP_PE_MAP_BAD_CHARACTER;
    // Check error_code
    switch (error_code) {
      case GT_IMP_PE_PENDING_BLOCKS:
      case GT_IMP_PE_MAP_PENDING_MAPS:
      case GT_IMP_PE_EOB:
      case GT_IMP_PE_MMAP_ATTRIBUTE_SCORE: break;
      default: // Sth went wrong
        if (map_parsed!=NULL) gt_map_delete(map_parsed);
        if (map_head!=NULL) gt_map_delete(map_head);
        return error_code;
    }
    // Check blocks & base lengths
    if (map_head==NULL) {
      map_head = map_parsed;
      last_map_block = map_head;
    } else {
      gt_map* const last_added = map_parsed;
      gt_map_set_next_block(last_map_block,last_added,QUIMERA,0);
      last_map_block = last_added;
    }
  } while (error_code==GT_IMP_PE_PENDING_BLOCKS);
  // Return Map
  *map = map_head;
  /*
   * Check error code
   */
  bool parse_attr = true, parsed_score = false;
  while (parse_attr) {
    switch (error_code) {
      case GT_IMP_PE_EOB:
      case GT_IMP_PE_MAP_PENDING_MAPS:
        parse_attr = false;
        break;
      case GT_IMP_PE_MMAP_ATTRIBUTE_SCORE:
        if (parsed_score) return GT_IMP_PE_MAP_BAD_CHARACTER;
        if ((error_code=gt_imp_parse_map_score_v1(text_line,&map_head->gt_score))) return error_code;
        map_head->phred_score = map_head->gt_score;
        parsed_score = true;
        error_code = gt_imp_next_element(text_line);
        break;
      case GT_IMP_PE_PENDING_BLOCKS: // No more pending blocks
        if (map_head!=NULL) gt_map_delete(map_head);
        return GT_IMP_PE_MAP_BAD_CHARACTER;
      default: // Sth went wrong
        if (map_head!=NULL) gt_map_delete(map_head);
        return error_code;
    }
  }
  return error_code;
}
GT_INLINE gt_status gt_imp_parse_alignment(
    const char** const text_line,gt_alignment* alignment,
    const bool has_quality_string,gt_map_parser_attributes* const map_parser_attr) {
  GT_NULL_CHECK(text_line);
  GT_ALIGNMENT_CHECK(alignment);
  gt_status error_code;
  // TAG
  if ((error_code=gt_imp_tag(text_line,alignment->tag,alignment->attributes))) return error_code;
  // READ
  error_code=gt_imp_read_block(text_line,alignment->read);
  if (gt_expect_false(error_code==GT_IMP_PE_PENDING_BLOCKS)) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
  if (gt_expect_false(error_code!=GT_IMP_PE_EOB)) return error_code;
  // QUALITIES
  if (has_quality_string) {
    if (gt_expect_false(**text_line==TAB)) {
      GT_NEXT_CHAR(text_line);
    } else {
      error_code=gt_imp_qualities_block(text_line,alignment->qualities);
      if (gt_expect_false(gt_string_get_length(alignment->qualities)!=
                          gt_string_get_length(alignment->read))) return GT_IMP_PE_QUAL_BAD_LENGTH;
      if (gt_expect_false(error_code==GT_IMP_PE_PENDING_BLOCKS)) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
      if (gt_expect_false(error_code!=GT_IMP_PE_EOB)) return error_code;
    }
  }
  // COUNTERS
  gt_vector_clear(gt_alignment_get_counters_vector(alignment));
  if ((error_code=gt_imp_counters(text_line,alignment->counters,alignment->attributes))) return error_code;
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  if (**text_line!=TAB) return GT_IMP_PE_BAD_SEPARATOR;
  GT_NEXT_CHAR(text_line);
  // MAPS
  error_code=gt_imp_parse_alignment_maps(text_line,alignment,map_parser_attr);
  return error_code;
}
GT_INLINE gt_status gt_imp_parse_template(
    const char** const text_line,gt_template* const template,
    const bool has_map_quality,gt_map_parser_attributes* const map_parser_attr) {
  GT_NULL_CHECK(text_line);
  GT_TEMPLATE_CHECK(template);
  gt_status error_code;
  // TAG
  if ((error_code=gt_imp_tag(text_line,template->tag,template->attributes))) {
    return error_code;
  }
  // READ
  uint64_t num_blocks = 0;
  error_code=GT_IMP_PE_PENDING_BLOCKS;
  while (error_code==GT_IMP_PE_PENDING_BLOCKS) {
    gt_alignment* const alignment = gt_template_get_block_dyn(template,num_blocks);
    error_code=gt_imp_read_block(text_line,alignment->read);
    if (error_code!=GT_IMP_PE_PENDING_BLOCKS && error_code!=GT_IMP_PE_EOB) return error_code;
    ++num_blocks;
  }
  // TAG Setup (Pair information based on template pair and num_blocks)
  gt_template_setup_pair_attributes_to_alignments(template,true);
  // QUALITIES
  if (has_map_quality) {
    if (gt_expect_false(**text_line==TAB)) {
      GT_NEXT_CHAR(text_line);
    } else {
      uint64_t i;
      error_code=GT_IMP_PE_PENDING_BLOCKS;
      for (i=0;i<num_blocks;++i) {
        if (error_code!=GT_IMP_PE_PENDING_BLOCKS) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
        gt_alignment* alignment = gt_template_get_block(template,i);
        error_code=gt_imp_qualities_block(text_line,alignment->qualities);
        if (gt_expect_false(gt_string_get_length(alignment->qualities)>0 &&
            gt_string_get_length(alignment->qualities)!=gt_string_get_length(alignment->read))){
          return GT_IMP_PE_QUAL_BAD_LENGTH;
        }
        if (error_code!=GT_IMP_PE_PENDING_BLOCKS && error_code!=GT_IMP_PE_EOB) return error_code;
      }
      if (error_code!=GT_IMP_PE_EOB) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
    }
  }
  // COUNTERS
  gt_vector_clear(gt_template_get_counters_vector(template));
  if (gt_expect_true(num_blocks>1)) {
    error_code=gt_imp_counters(text_line,gt_template_get_counters_vector(template),template->attributes);
  } else {
    gt_alignment* const alignment = gt_template_get_block(template,0);
    error_code=gt_imp_counters(text_line,gt_alignment_get_counters_vector(alignment),alignment->attributes);
  }
  if (error_code) return error_code;
  if (gt_expect_false((**text_line)!=TAB)) return GT_IMP_PE_PREMATURE_EOL;
  GT_NEXT_CHAR(text_line);
  // MAPS
  if (gt_expect_true(num_blocks>1)) {
    error_code = gt_imp_parse_template_maps(text_line,template,map_parser_attr);
  } else {
    error_code = gt_imp_parse_alignment_maps(text_line,gt_template_get_block(template,0),map_parser_attr);
  }
  return error_code;
}
/*
 * MAP string parsers
 */
GT_INLINE gt_status gt_input_map_parse_counters(const char* const string,gt_vector* const counters,gt_attributes* attributes) {
  GT_NULL_CHECK(string);
  GT_VECTOR_CHECK(counters);
  GT_ATTRIBUTES_CHECK(attributes);
  const char* _string = string; // Placeholder
  return gt_imp_counters(&_string,counters,attributes);
}
GT_INLINE gt_status gt_input_map_parse_map(const char* const string,gt_map** const map,gt_map_parser_attributes* map_parser_attr) {
  GT_NULL_CHECK(string);
  GT_NULL_CHECK(map);
  GT_MAP_PARSER_CHECK_ATTRIBUTES(map_parser_attr);
  // Placeholders
  const char* _string = string;
  gt_status error_code = gt_imp_map_blocks(&_string,map,map_parser_attr);
  GT_IMP_PARSE_MAP_BEGIN_ERROR_HANDLER(error_code) {
    return error_code;
  } GT_IMP_PARSE_MAP_END_ERROR_HANDLER(error_code)
  return 0;
}
GT_INLINE gt_status gt_input_map_parse_map_list(const char* const string,gt_vector* const maps,gt_map_parser_attributes* map_parser_attr) {
  GT_NULL_CHECK(string);
  GT_VECTOR_CHECK(maps);
  GT_MAP_PARSER_CHECK_ATTRIBUTES(map_parser_attr);
  const char* _string = string; // Placeholder
  // Read maps
  uint64_t num_parsed_maps = 0;
  gt_map* map = NULL;
  gt_status error_code=GT_IMP_PE_MAP_PENDING_MAPS;
  while (error_code==GT_IMP_PE_MAP_PENDING_MAPS && num_parsed_maps<map_parser_attr->max_parsed_maps) {
    // Parse map list
    gt_status error_code = gt_imp_map_blocks(&_string,&map,map_parser_attr);
    switch (error_code) {
      case GT_IMP_PE_MAP_PENDING_MAPS:
      case GT_IMP_PE_EOB:
        break;
      default:
        return error_code;
        break;
    }
    // Add map (NO check duplicates)
    gt_vector_insert(maps,map,gt_map*);
    ++num_parsed_maps;
  }
  GT_IMP_PARSE_MAP_BEGIN_ERROR_HANDLER(error_code) {
    return error_code;
  } GT_IMP_PARSE_MAP_END_ERROR_HANDLER(error_code)
  return 0;
}
GT_INLINE gt_status gt_input_map_parse_alignment(const char* const string,gt_alignment* const alignment) {
  GT_NULL_CHECK(string);
  GT_ALIGNMENT_CHECK(alignment);
  const char* _string = string; // Placeholder
  gt_alignment_clear(alignment); // Clear alignment
  // Count fields
  uint64_t num_fields=0, i=0;
  while (gt_expect_true(_string[i]!=EOS)) {
    if (gt_expect_false(_string[i]==TAB)) ++num_fields;
    i++;
  }
  if (gt_expect_true(num_fields==4 || num_fields==3)) {
    gt_map_parser_attributes map_parser_attr = GT_MAP_PARSER_ATTR_DEFAULT(false);
    return gt_imp_parse_alignment(&_string,alignment,num_fields==4,&map_parser_attr);
  }
  return GT_IMP_PE_WRONG_FILE_FORMAT;
}
GT_INLINE gt_status gt_input_map_parse_template(const char* const string,gt_template* const template) {
  GT_NULL_CHECK(string);
  GT_TEMPLATE_CHECK(template);
  const char* _string = string; // Placeholder
  gt_template_clear(template,true); // Clear template
  // Count fields
  uint64_t num_fields=0, i=0;
  while (gt_expect_true(_string[i]!=EOS)) {
    if (gt_expect_false(_string[i]==TAB)) ++num_fields;
    i++;
  }
  if (gt_expect_true(num_fields==4 || num_fields==3)) {
    gt_map_parser_attributes map_parser_attr = GT_MAP_PARSER_ATTR_DEFAULT(false);
    return gt_imp_parse_template(&_string,template,num_fields==4,&map_parser_attr);
  }
  return GT_IMP_PE_WRONG_FILE_FORMAT;
}
/*
 * MAP High-level Parsers
 */
GT_INLINE gt_status gt_imp_get_template(
    gt_buffered_input_file* const buffered_map_input,
    gt_template* const template,gt_map_parser_attributes* const map_parser_attr) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(map_parser_attr);
  gt_input_file* const input_file = buffered_map_input->input_file;
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_map_input)) {
    if ((error_code=gt_input_map_parser_reload_buffer(buffered_map_input,true,GT_IMP_NUM_LINES))!=GT_INPUT_STATUS_OK) return error_code;
  }
  // Check file format
  if (gt_input_map_parser_check_map_file_format(buffered_map_input)) {
    gt_fatal_error(PARSE_MAP_BAD_FILE_FORMAT,gt_input_file_get_file_name(input_file),buffered_map_input->current_line_num);
    return GT_INPUT_STATUS_FAIL;
  }
  // Prepare the template
  char* const line_start = buffered_map_input->cursor;
  const uint64_t line_num = buffered_map_input->current_line_num;
  gt_template_clear(template,true);
  template->template_id = line_num;
  // Parse template
  if ((error_code=gt_imp_parse_template((const char** const)&(buffered_map_input->cursor),
      template,input_file->map_type.contains_qualities,map_parser_attr))) {
    gt_input_map_parser_prompt_error(buffered_map_input,line_num,
        buffered_map_input->cursor-line_start,error_code);
    gt_input_map_parser_next_record(buffered_map_input);
    if (map_parser_attr->src_text!=NULL) {
      gt_string_set_nstring(map_parser_attr->src_text,line_start,buffered_map_input->cursor-line_start);
    }
    return GT_INPUT_STATUS_FAIL;
  }
  // Store source record
  if (map_parser_attr->src_text!=NULL) {
    gt_string_set_nstring(map_parser_attr->src_text,line_start,buffered_map_input->cursor-line_start);
  }
  // Next record
  gt_input_map_parser_next_record(buffered_map_input);
  return GT_INPUT_STATUS_OK;
}
GT_INLINE gt_status gt_imp_get_alignment(
    gt_buffered_input_file* const buffered_map_input,
    gt_alignment* const alignment,gt_map_parser_attributes* const map_parser_attr) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(map_parser_attr);
  gt_input_file* const input_file = buffered_map_input->input_file;
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_map_input)) {
    if ((error_code=gt_input_map_parser_reload_buffer(buffered_map_input,false,GT_IMP_NUM_LINES))!=GT_INPUT_STATUS_OK) return error_code;
  }
  // Check file format
  if (gt_input_map_parser_check_map_file_format(buffered_map_input)) {
    gt_error(PARSE_MAP_BAD_FILE_FORMAT,gt_input_file_get_file_name(input_file),buffered_map_input->current_line_num);
    return GT_INPUT_STATUS_FAIL;
  }
  // Allocate memory for the alignment
  char* const line_start = buffered_map_input->cursor;
  const uint64_t line_num = buffered_map_input->current_line_num;
  gt_alignment_clear(alignment);
  alignment->alignment_id = line_num;
  // Parse alignment
  if ((error_code=gt_imp_parse_alignment((const char** const)&(buffered_map_input->cursor),
      alignment,input_file->map_type.contains_qualities,map_parser_attr))) {
    gt_input_map_parser_prompt_error(buffered_map_input,line_num,
        buffered_map_input->cursor-line_start,error_code);
    gt_input_map_parser_next_record(buffered_map_input);
    if (map_parser_attr->src_text) {
      gt_string_set_nstring(map_parser_attr->src_text,line_start,buffered_map_input->cursor-line_start);
    }
    return GT_INPUT_STATUS_FAIL;
  }
  // Store source record
  if (map_parser_attr->src_text) {
    gt_string_set_nstring(map_parser_attr->src_text,line_start,buffered_map_input->cursor-line_start);
  }
  // Next record
  gt_input_map_parser_next_record(buffered_map_input);
  return GT_INPUT_STATUS_OK;
}
/*
 * MAP High-level Parsers
 *   - High-level parsing to extract one template/alignment from the buffered file (reads one line)
 *   - Syntax checking
 *   - Transparent buffer block reload
 *   - Template/Alignment transparent memory management
 */
GT_INLINE gt_status gt_input_map_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,gt_map_parser_attributes* map_parser_attr) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  GT_MAP_PARSER_CHECK_ATTRIBUTES(map_parser_attr);
  gt_status error_code;
  if ((error_code=gt_imp_get_template(buffered_map_input,template,map_parser_attr))!=GT_INPUT_STATUS_OK) {
    return (error_code==GT_INPUT_STATUS_EOF) ? GT_INPUT_STATUS_EOF : GT_INPUT_STATUS_FAIL;
  }
  if (gt_template_get_num_blocks(template)==1 && map_parser_attr->force_read_paired) {
    if ((error_code=gt_imp_get_alignment(buffered_map_input,gt_template_get_block_dyn(template,1),map_parser_attr))!=GT_INPUT_STATUS_OK) {
      return GT_INPUT_STATUS_FAIL;
    }
    // Check TAG consistency
    gt_alignment* const end1 = gt_template_get_block(template,0);
    gt_alignment* const end2 = gt_template_get_block(template,1);
    if (!gt_string_equals(end1->tag,end2->tag)) return GT_INPUT_STATUS_FAIL;
    // TAG Setup
    gt_template_setup_pair_attributes_to_alignments(template,false);
  }
  return error_code;
}
GT_INLINE gt_status gt_input_map_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,gt_map_parser_attributes* map_parser_attr) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  GT_MAP_PARSER_CHECK_ATTRIBUTES(map_parser_attr);
  return gt_imp_get_alignment(buffered_map_input,alignment,map_parser_attr);
}
/*
 * Synch read of blocks
 */
GT_INLINE gt_status gt_input_map_parser_synch_blocks(
    gt_buffered_input_file* const buffered_input1,gt_buffered_input_file* const buffered_input2,pthread_mutex_t* const input_mutex) {
  gt_map_parser_attributes map_parser_attr = GT_MAP_PARSER_ATTR_DEFAULT(true);
  return gt_input_map_parser_synch_blocks_va(input_mutex,&map_parser_attr,2,buffered_input1,buffered_input2);
}
GT_INLINE gt_status gt_input_map_parser_synch_blocks_v(
    pthread_mutex_t* const input_mutex,gt_map_parser_attributes* const map_parser_attr,
    uint64_t num_inputs,gt_buffered_input_file* const buffered_input,va_list v_args) {
  GT_NULL_CHECK(input_mutex);
  GT_NULL_CHECK(map_parser_attr);
  GT_ZERO_CHECK(num_inputs);
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed (synch)
  if (gt_buffered_input_file_eob(buffered_input)) {
    // Dump buffer if BOF it attached to the input, and get new out block (always FIRST)
    gt_buffered_input_file_dump_attached_buffers(buffered_input->attached_buffered_output_file);
    // Read synch blocks & Reload all the 'buffered_input' files
    GT_BEGIN_MUTEX_SECTION(*input_mutex) {
      if ((error_code=gt_input_map_parser_reload_buffer(buffered_input,
          map_parser_attr->force_read_paired,GT_IMP_NUM_LINES))!=GT_INPUT_STATUS_OK) {
        GT_END_MUTEX_SECTION(*input_mutex);
        return error_code;
      }
      // Reload the rest of the 'buffered_input' files
      while ((--num_inputs)>0) {
        gt_buffered_input_file* remaining_buffered_input = va_arg(v_args,gt_buffered_input_file*);
        GT_BUFFERED_INPUT_FILE_CHECK(remaining_buffered_input);
        if ((error_code=gt_input_map_parser_reload_buffer(remaining_buffered_input,
            map_parser_attr->force_read_paired,GT_IMP_NUM_LINES))!=GT_INPUT_STATUS_OK) {
          GT_END_MUTEX_SECTION(*input_mutex);
          return error_code;
        }
      }
    } GT_END_MUTEX_SECTION(*input_mutex);
  }
  return GT_INPUT_STATUS_OK;
}
GT_INLINE gt_status gt_input_map_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,gt_map_parser_attributes* const map_parser_attr,
    const uint64_t num_inputs,gt_buffered_input_file* const buffered_input,...) {
  GT_NULL_CHECK(input_mutex);
  GT_NULL_CHECK(map_parser_attr);
  GT_ZERO_CHECK(num_inputs);
  va_list v_args;
  va_start(v_args,buffered_input);
  gt_status error_code =
      gt_input_map_parser_synch_blocks_v(input_mutex,map_parser_attr,num_inputs,buffered_input,v_args);
  va_end(v_args);
  return error_code;
}
GT_INLINE gt_status gt_input_map_parser_synch_blocks_a(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file** const buffered_input,
    const uint64_t num_inputs,gt_map_parser_attributes* const map_parser_attr) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input[0]);
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed (synch)
  if (gt_buffered_input_file_eob(buffered_input[0])) {
    /*
     * Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
     */
    gt_buffered_input_file_dump_attached_buffers(buffered_input[0]->attached_buffered_output_file);
    /*
     * Read synch blocks
     */
    GT_BEGIN_MUTEX_SECTION(*input_mutex) {
      uint64_t i;
      for (i=0;i<num_inputs;++i) {
        // Reload the 'buffered_input' files
        GT_BUFFERED_INPUT_FILE_CHECK(buffered_input[i]);
        if ((error_code=gt_input_map_parser_reload_buffer(buffered_input[i],
            map_parser_attr->force_read_paired,GT_IMP_NUM_LINES))!=GT_INPUT_STATUS_OK) {
          GT_END_MUTEX_SECTION(*input_mutex);
          return error_code;
        }
      }
    } GT_END_MUTEX_SECTION(*input_mutex);
  }
  return GT_INPUT_STATUS_OK;
}
// Used to merge files in parallel
GT_INLINE gt_status gt_input_map_parser_synch_blocks_by_subset(
    pthread_mutex_t* const input_mutex,gt_map_parser_attributes* const map_parser_attr,
    gt_buffered_input_file* const buffered_map_input_master,gt_buffered_input_file* const buffered_map_input_slave) {
  GT_NULL_CHECK(input_mutex);
  GT_NULL_CHECK(map_parser_attr);
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input_master);
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input_slave);
  gt_status error_code_master, error_code_slave;
  // Check the end_of_block. Reload buffer if needed (synch)
  const bool eob_master = gt_buffered_input_file_eob(buffered_map_input_master);
  const bool eob_slave = gt_buffered_input_file_eob(buffered_map_input_slave);
  if (!eob_master) return GT_INPUT_STATUS_OK;
  if (!eob_slave) return GT_INPUT_STATUS_FAIL;
  /*
   * Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
   */
  gt_buffered_input_file_dump_attached_buffers(buffered_map_input_master->attached_buffered_output_file);
  /*
   * Read synch blocks
   */
  GT_BEGIN_MUTEX_SECTION(*input_mutex) {
    /*
     * Read new input block for the slave
     */
    if ((error_code_slave=gt_input_map_parser_reload_buffer(buffered_map_input_slave,
        map_parser_attr->force_read_paired,GT_IMP_SUBSET_NUM_LINES+1))!=GT_INPUT_STATUS_OK) {
      // Read for the master and finish
      error_code_master=gt_input_map_parser_reload_buffer(buffered_map_input_master,
          map_parser_attr->force_read_paired,GT_IMP_SUBSET_NUM_LINES+1);
      GT_END_MUTEX_SECTION(*input_mutex);
      return error_code_master;
    }
    /*
     * Read new input block for the master
     */
    // Get the last tag of the slave block
    gt_string* const last_tag = gt_string_new(0);
    if ((error_code_slave=gt_imp_get_tag_last_read(buffered_map_input_slave,last_tag)!=GT_INPUT_STATUS_OK)) {
      GT_END_MUTEX_SECTION(*input_mutex);
      return error_code_slave;
    }
    // Read new input block for master matching the last tag of the slave block
    if ((error_code_master=gt_imp_reload_buffer_matching_tag(
        buffered_map_input_master,last_tag,map_parser_attr->force_read_paired))!=GT_INPUT_STATUS_OK) {
      gt_string_delete(last_tag); // Free
      GT_END_MUTEX_SECTION(*input_mutex);
      return error_code_master;
    }
    gt_string_delete(last_tag); // Free
  } GT_END_MUTEX_SECTION(*input_mutex);
  return GT_INPUT_STATUS_OK;
}


