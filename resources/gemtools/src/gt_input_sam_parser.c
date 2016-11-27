/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_sam_parser.c
 * DATE: 17/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for SAM format
 */

#include "gt_input_sam_parser.h"
#include "gt_input_parser.h"

// Constants
#define GT_ISP_NUM_LINES GT_NUM_LINES_10K
#define GT_ISP_NUM_INITIAL_MAPS 5

/*
 * SAM parser attributes
 */
GT_INLINE gt_sam_parser_attributes* gt_input_sam_parser_attributes_new() {
  gt_sam_parser_attributes* attributes = gt_alloc(gt_sam_parser_attributes);
  gt_input_sam_parser_attributes_reset_defaults(attributes);
  return attributes;
}
GT_INLINE void gt_input_sam_parser_attributes_delete(gt_sam_parser_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  free(attributes);
}
GT_INLINE void gt_input_sam_parser_attributes_reset_defaults(gt_sam_parser_attributes* const attributes) {
  attributes->sam_soap_style = false;
}
GT_INLINE void gt_input_sam_parser_attributes_set_soap_compilant(gt_sam_parser_attributes* const attributes) {
  attributes->sam_soap_style = true;
}

// Internal pair-pending
typedef struct {
  // Current map info
  gt_string map_seq_name;
  uint64_t map_position;
  uint64_t end_position; // 0/1
  // Next map info
  gt_string next_seq_name;
  uint64_t next_position;
  // Map location and span info
  uint64_t map_displacement; // In alignment's map vector
  uint64_t num_maps; // Maps in the vector coupled to the first one
} gt_sam_pending_end;

#define GT_SAM_INIT_PENDING { .map_seq_name.allocated=0, .next_seq_name.allocated=0 }

/*
 * SAM File Format test
 */
#define GT_INPUT_FILE_SAM_READ_HEADERS_CMP_TAG(tag_array,l1,l2) ((tag_array)[0]==l1 && (tag_array)[1]==l2 && (tag_array)[2]==TAB)
#define GT_INPUT_FILE_SAM_READ_HEADERS_CMP_ATTR(tag_array,l1,l2) ((tag_array)[0]==l1 && (tag_array)[1]==l2 && (tag_array)[2]==COLON)
GT_INLINE gt_status gt_input_file_sam_read_headers(
    char* const buffer,const uint64_t buffer_size,gt_sam_headers* const sam_headers,
    uint64_t* const characters_read,uint64_t* const lines_read) {
  uint64_t buffer_pos=0, lines=0;
  // Read until no more header lines are parsed
  while (buffer[buffer_pos]==GT_SAM_HEADER_BEGIN) {
    ++buffer_pos;
    if (GT_INPUT_FILE_SAM_READ_HEADERS_CMP_TAG(buffer+buffer_pos,'H','D')) {
      buffer_pos+=3;
      while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL) ++buffer_pos;
      if (buffer_pos==buffer_size) return 0;
      ++buffer_pos;
    } else if (GT_INPUT_FILE_SAM_READ_HEADERS_CMP_TAG(buffer+buffer_pos,'S','Q')) {
      buffer_pos+=3;
      while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL) ++buffer_pos;
      if (buffer_pos==buffer_size) return 0;
      ++buffer_pos;
    } else if (GT_INPUT_FILE_SAM_READ_HEADERS_CMP_TAG(buffer+buffer_pos,'R','G')) {
      buffer_pos+=3;
      while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL) ++buffer_pos;
      if (buffer_pos==buffer_size) return 0;
      ++buffer_pos;
    } else if (GT_INPUT_FILE_SAM_READ_HEADERS_CMP_TAG(buffer+buffer_pos,'P','G')) {
      buffer_pos+=3;
      while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL) ++buffer_pos;
      if (buffer_pos==buffer_size) return 0;
      ++buffer_pos;
    } else if (GT_INPUT_FILE_SAM_READ_HEADERS_CMP_TAG(buffer+buffer_pos,'C','O')) {
      buffer_pos+=3;
      while (buffer_pos<buffer_size && buffer[buffer_pos]!=EOL) ++buffer_pos;
      if (buffer_pos==buffer_size) return 0;
      ++buffer_pos;
    } else {
      return -1;
    }
    ++lines;
  }
  *characters_read = buffer_pos;
  *lines_read = lines;
  return 0;
}

#define GT_ISP_TEST_SAM_SKIP_STRING() \
  while (buffer_pos<buffer_size && buffer[buffer_pos]!=TAB && buffer[buffer_pos]!=EOL) ++buffer_pos; \
  if (buffer_pos==buffer_size || buffer[buffer_pos]==EOL) return false; \
  ++buffer_pos

#define GT_ISP_TEST_SAM_TEST_INTEGER() \
  while (buffer_pos<buffer_size && buffer[buffer_pos]!=TAB && buffer[buffer_pos]!=EOL) { \
    if (!gt_is_number(buffer[buffer_pos])) return false; \
    ++buffer_pos; \
  } \
  if (buffer_pos==buffer_size || buffer[buffer_pos]==EOL) return false; \
  ++buffer_pos

GT_INLINE bool gt_input_sam_parser_test_sam(
    char* const file_name,const uint64_t line_num,char* const buffer,const uint64_t buffer_size,
    uint64_t* const characters_read,uint64_t* const lines_read,gt_sam_headers* const sam_headers,const bool show_errors) {
  /*
   * (1) @SQ     SN:chr10        LN:135534747
   *     @SQ     SN:chr11        LN:135006516
   * (2) 1/1     16  chr12  57338496  37  75M  *  0  0  TCTG...TTGN  BbQ..QQB  XT:A:U  NM:i:1
   */
  uint64_t buffer_pos=0;
  if (buffer[0]==GT_SAM_HEADER_BEGIN) { // Read headers
    if(gt_input_file_sam_read_headers(buffer,buffer_size,sam_headers,characters_read,lines_read)!=0) return false;
    buffer_pos = *characters_read;
  } else {
    *lines_read = 0;
    *characters_read = 0;
  }
  /*
   * Check SAM record
   *   SRR003161.1     4       *       0       0       *       *       0       0       *       *       AS:i:0
   */
  // Skip TAG
  GT_ISP_TEST_SAM_SKIP_STRING();
  // Check FLAG
  GT_ISP_TEST_SAM_TEST_INTEGER();
  // Skip RNAME
  GT_ISP_TEST_SAM_SKIP_STRING();
  // Check POS
  GT_ISP_TEST_SAM_TEST_INTEGER();
  // Check MAPQ
  GT_ISP_TEST_SAM_TEST_INTEGER();
  // Skip CIGAR
  GT_ISP_TEST_SAM_SKIP_STRING();
  // Skip RNEXT
  GT_ISP_TEST_SAM_SKIP_STRING();
  // Check PNEXT
  GT_ISP_TEST_SAM_TEST_INTEGER();
  // Check TLEN
  GT_ISP_TEST_SAM_SKIP_STRING();
  // Check SEQ
  GT_ISP_TEST_SAM_SKIP_STRING();
  // Skip QUAL
  while (buffer_pos<buffer_size && buffer[buffer_pos]!=TAB && buffer[buffer_pos]!=EOL) ++buffer_pos;
  if (buffer_pos==buffer_size) return false;
  return true;
}

#define GT_ISP_HEADERS_WRONG_FORMAT -1
#define GT_ISP_HEADERS_END 0

GT_INLINE bool gt_input_file_test_sam(
    gt_input_file* const input_file,gt_sam_headers* const sam_headers,const bool show_errors) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_NULL_CHECK(sam_headers);
  uint64_t characters_read = 0, processed_lines = 0;
  if (gt_input_sam_parser_test_sam(
      gt_input_file_get_file_name(input_file),input_file->processed_lines+1,
      (char*)input_file->file_buffer,input_file->buffer_size,&characters_read,&processed_lines,sam_headers,show_errors)) {
    input_file->buffer_begin = characters_read;
    input_file->buffer_pos = characters_read;
    input_file->processed_lines = processed_lines;
    return true;
  }
  return false;
}
GT_INLINE gt_status gt_input_sam_parser_check_sam_file_format(gt_buffered_input_file* const buffered_sam_input) {
  gt_input_file* const input_file = buffered_sam_input->input_file;
  if (gt_expect_false(input_file->file_format==FILE_FORMAT_UNKNOWN)) { // Unknown
    gt_sam_headers sam_headers;
    // Mutex format detection (because the first one must read the headers)
    gt_input_file_lock(input_file);
      const bool is_sam_format =
          gt_input_file_test_sam(input_file,&sam_headers,true);
    gt_input_file_unlock(input_file);
    if (!is_sam_format) return GT_ISP_PE_WRONG_FILE_FORMAT;
    input_file->file_format = SAM;
  } else if (gt_expect_false(input_file->file_format!=SAM)) {
    return GT_ISP_PE_WRONG_FILE_FORMAT;
  }
  return 0;
}

/*
 * SAM File basics
 */
/* Error handler */
GT_INLINE void gt_input_sam_parser_prompt_error(
    gt_buffered_input_file* const buffered_sam_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code) {
  // Display textual error msg
  const char* const file_name = (buffered_sam_input != NULL) ?
      gt_input_file_get_file_name(buffered_sam_input->input_file) : "<<LazyParsing>>";
  if ((buffered_sam_input == NULL)) {
    line_num = 0; column_pos = 0;
  }
  switch (error_code) {
    case 0: /* No error */ break;
    case GT_ISP_PE_WRONG_FILE_FORMAT: gt_error(PARSE_SAM_BAD_FILE_FORMAT,file_name,line_num,column_pos); break;
    case GT_ISP_PE_PREMATURE_EOL: gt_error(PARSE_SAM_PREMATURE_EOL,file_name,line_num,column_pos); break;
    case GT_ISP_PE_EXPECTED_NUMBER: gt_error(PARSE_SAM_EXPECTED_NUMBER,file_name,line_num,column_pos); break;
    case GT_ISP_PE_BAD_CHARACTER: gt_error(PARSE_SAM_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_ISP_PE_WRONG_READ_CONTENT: gt_error(PARSE_SAM_WRONG_READ_CONTENT,file_name,line_num,column_pos); break;
    case GT_ISP_PE_CIGAR_PREMATURE_END: gt_error(PARSE_SAM_CIGAR_PREMATURE_END,file_name,line_num,column_pos); break;
    case GT_ISP_PE_SAM_UNMAPPED_XA: gt_error(PARSE_SAM_UNMAPPED_XA,file_name,line_num,column_pos); break;
    case GT_ISP_PE_WRONG_NUM_XA: gt_error(PARSE_SAM_WRONG_NUM_XA,file_name,line_num,column_pos); break;
    case GT_ISP_PE_UNSOLVED_PENDING_MAPS: gt_error(PARSE_SAM_UNSOLVED_PENDING_MAPS,file_name,line_num,column_pos); break;
    default:
      gt_error(PARSE_SAM,file_name,line_num,column_pos);
      break;
  }
}
/* MAP file. Skip record */
GT_INLINE void gt_input_sam_parser_next_record(gt_buffered_input_file* const buffered_sam_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  if (!gt_buffered_input_file_eob(buffered_sam_input)) {
    gt_buffered_input_file_skip_line(buffered_sam_input);
  }
}
/*
 * SAM file. Reload internal buffer
 */
/* SAM file. Synchronized get block wrt to sam records */
GT_INLINE gt_status gt_input_sam_parser_get_block(
    gt_buffered_input_file* const buffered_sam_input,const uint64_t num_records) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  gt_input_file* const input_file = buffered_sam_input->input_file;
  // Read lines
  if (gt_input_file_eof(input_file)) return GT_INPUT_STATUS_EOF;
  gt_input_file_lock(input_file);
  if (gt_input_file_eof(input_file)) {
    gt_input_file_unlock(input_file);
    return GT_INPUT_STATUS_EOF;
  }
  buffered_sam_input->block_id = gt_input_file_get_next_id(input_file);
  buffered_sam_input->current_line_num = input_file->processed_lines+1;
  gt_vector_clear(buffered_sam_input->block_buffer); // Clear dst buffer
  // Read lines & synch SAM records
  uint64_t lines_read = 0;
  while (lines_read<num_records &&
      gt_input_file_next_line(input_file,buffered_sam_input->block_buffer) ) ++lines_read;
  if (lines_read==num_records) { // !EOF, Synch wrt to tag content
    uint64_t num_blocks=0, num_tabs=0;
    gt_string* const reference_tag = gt_string_new(30);
    if (gt_input_file_next_record(input_file,buffered_sam_input->block_buffer,reference_tag,&num_blocks,&num_tabs)) {
      gt_input_parse_tag_chomp_pairend_info(reference_tag);
      while (gt_input_file_next_record_cmp_first_field(input_file,reference_tag)) {
        if (!gt_input_file_next_record(input_file,buffered_sam_input->block_buffer,NULL,&num_blocks,&num_tabs)) break;
        ++lines_read;
      }
    }
    gt_string_delete(reference_tag);
  }
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffered_sam_input->block_buffer);
  if (lines_read > 0 && *gt_vector_get_last_elm(buffered_sam_input->block_buffer,char) != EOL) {
    gt_vector_insert(buffered_sam_input->block_buffer,EOL,char);
  }
  input_file->processed_lines+=lines_read;
  buffered_sam_input->lines_in_buffer = lines_read;
  gt_input_file_unlock(input_file);

  // Setup the block
  buffered_sam_input->cursor = gt_vector_get_mem(buffered_sam_input->block_buffer,char);
  return buffered_sam_input->lines_in_buffer;
}
/* SAM file. Reload internal buffer */
GT_INLINE gt_status gt_input_sam_parser_reload_buffer(gt_buffered_input_file* const buffered_sam_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  // Dump buffer if BOF it attached to SAM-input, and get new out block (always FIRST)
  gt_buffered_input_file_dump_attached_buffers(buffered_sam_input->attached_buffered_output_file);
  // Read new input block
  const uint64_t read_lines =
      gt_input_sam_parser_get_block(buffered_sam_input,GT_ISP_NUM_LINES);
  if (gt_expect_false(read_lines==0)) return GT_ISP_EOF;
  // Assign block ID
  gt_buffered_input_file_set_id_attached_buffers(buffered_sam_input->attached_buffered_output_file,buffered_sam_input->block_id);
  return GT_ISP_OK;
}
/*
 * SAM format. Basic building block for parsing
 */
GT_INLINE uint64_t gt_isp_read_tag(const char** const init_text_line,const char** const end_text_line,gt_string* const tag) {
  // Save text_line state
  const char* text_cp = *init_text_line;
  const char** const ptext_cp = &text_cp;
  // Read tag
  const char* const tag_begin = *ptext_cp;
  GT_READ_UNTIL(ptext_cp,**ptext_cp==TAB || **ptext_cp==SPACE);
  if (GT_IS_EOL(ptext_cp)) return GT_ISP_PE_PREMATURE_EOL;
  // Set tag
  uint64_t const tag_length = *ptext_cp-tag_begin;
  gt_string_set_nstring_static(tag,tag_begin,tag_length);
  // Read the rest till next field
  if (**ptext_cp==SPACE) {
    GT_READ_UNTIL(ptext_cp,**ptext_cp==TAB);
    if (GT_IS_EOL(ptext_cp)) return GT_ISP_PE_PREMATURE_EOL;
  }
  GT_NEXT_CHAR(ptext_cp);
  *end_text_line = *ptext_cp;
  return 0;
}
/*
 * SAM CIGAR ::
 *   2M503N34M757N40M || 5M1D95M3I40M || ...
 */
GT_INLINE gt_status gt_isp_parse_sam_cigar(char** const text_line,gt_map** _map,const bool reverse_strand) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_NULL_CHECK(_map); GT_MAP_CHECK(*_map);
  gt_map* map = *_map;
  // Clear mismatches
  gt_map_clear_misms(map);
  if (**text_line==STAR) { // No CIGAR available
    GT_NEXT_CHAR(text_line);
    return 0;
  }
  // Aux variables as to track the position in the read and the genome span
  uint64_t length, position = 0, reference_span=0;
  while (**text_line!=TAB && **text_line!=EOL) {
    // Parse misms_op length
    if (!gt_is_number(**text_line)) return GT_ISP_PE_EXPECTED_NUMBER;
    GT_PARSE_NUMBER(text_line,length);
    // Parse misms_op
    if (gt_expect_false(**text_line==EOL || **text_line==TAB)) return GT_ISP_PE_CIGAR_PREMATURE_END;
    gt_misms misms;
    const char cigar_op = **text_line;
    GT_NEXT_CHAR(text_line);
    switch (cigar_op) {
      case 'M':
      case '=':
      case 'X':
        position += length;
        reference_span += length;
        break;
      case 'P': // Padding. Nothing specific implemented
      case 'S': // Soft clipping. Nothing specific implemented
      case 'H': // Hard clipping. Nothing specific implemented (we don't even store this)
        // break; //FIXME
      case 'I': // Insertion to the reference
        misms.misms_type = DEL;
        misms.position = position;
        misms.size = length;
        position += length;
        gt_map_add_misms(map,&misms);
        break;
      case 'D': // Deletion from the reference
        misms.misms_type = INS;
        misms.position = position;
        misms.size = length;
        reference_span += length;
        gt_map_add_misms(map,&misms);
        break;
      case 'N': { // Split. Eg TOPHAT, GEM, ...
        // Create a new map block
        gt_map* next_map = gt_map_new();
        gt_map_set_seq_name(next_map,gt_map_get_seq_name(map),gt_map_get_seq_name_length(map));
        gt_map_set_position(next_map,gt_map_get_position(map)+reference_span+length);
        gt_map_set_strand(next_map,gt_map_get_strand(map));
        gt_map_set_base_length(next_map,gt_map_get_base_length(map)-position);
        // Close current map block
        gt_map_set_base_length(map,position);
        if (reverse_strand) {
          gt_map_set_next_block(next_map,map,SPLICE,length);
        } else {
          gt_map_set_next_block(map,next_map,SPLICE,length);
        }
        // Swap maps & Reset position,reference_span
        map = next_map;
        position=0; reference_span=0;
        }
        break;
      default:
        return GT_ISP_PE_BAD_CHARACTER;
        break;
    }
  }
  gt_map_set_base_length(map,position);
//  // Consider map CIGAR in the reverse strand
//  if (reverse_strand) {
//    *_map = map;
//    GT_MAP_ITERATE(map,map_it) {
//      gt_map_reverse_misms(map_it);
//    }
//  }
  return 0;
}

#define GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL() \
  if (GT_IS_EOL(text_line)) { \
    gt_map_delete(map); return GT_ISP_PE_PREMATURE_EOL; \
  }
#define GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT() \
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL(); \
  GT_NEXT_CHAR(text_line)
#define GT_ISP_PARSE_SAM_ALG_SKIP_FIELD() \
  GT_READ_UNTIL(text_line,**text_line==TAB); \
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT()
#define GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(number) \
  if (!gt_is_number(**text_line)) { gt_map_delete(map); return GT_ISP_PE_EXPECTED_NUMBER; } \
  GT_PARSE_NUMBER(text_line,number)
#define GT_ISP_IF_OPT_FIELD(text_line,char1,char2,type_char) { \
  if ((*text_line)[0]==char1 && (*text_line)[1]==char2 && (*text_line)[2]!=EOL && (*text_line)[3]==type_char) {
#define GT_ISP_END_OPT_FIELD }}

GT_INLINE void gt_isp_parse_sam_opt_md_next_mismatch(
    gt_map* const map,
    gt_misms** const misms,
    uint64_t* const position) {
  if (*misms==NULL) {
    *position = 0;
  } else {
    ++(*position);
  }
  if (*position < gt_map_get_num_misms(map)) {
    *misms = gt_map_get_misms(map,*position);
  } else {
    *position = UINT64_MAX;
    *misms = NULL;
  }
}
GT_INLINE gt_status gt_isp_parse_sam_opt_md(
    const char** const text_line,
    gt_alignment* const alignment,
    gt_vector* const maps_vector,
    gt_sam_pending_end* const pending) {
  // Skip TAG
  *text_line+=5;
  // Need to begin with number
  if (!gt_is_number(**text_line)) return GT_ISP_PE_EXPECTED_NUMBER;
  // Allocate new misms vector
  gt_vector* misms_vector = gt_vector_new(4,sizeof(gt_misms));
  // Prepare Mismatches-iterator
  gt_map* const map = *gt_vector_get_elm(maps_vector,0,gt_map*);
  gt_misms* misms = NULL;
  uint64_t misms_pos;
  gt_isp_parse_sam_opt_md_next_mismatch(map,&misms,&misms_pos);
  // Parse Mismatches
  uint64_t key_pos=0, text_pos=0, number=0;
  while (**text_line!=TAB && **text_line!=EOL) {
    // Check MD CIGAR
    if (**text_line=='^') {
      ++(*text_line);
      while (gt_is_dna(**text_line)) ++(*text_line);
    } else if (gt_is_dna(**text_line)) {
      gt_misms mismatch;
      mismatch.misms_type = MISMS;
      mismatch.position = key_pos;
      mismatch.base = **text_line;
      gt_vector_insert(misms_vector,mismatch,gt_misms);
      ++(*text_line);
      ++key_pos;
    } else {
      // Match
      if (!gt_is_number(**text_line)) {
        fprintf(stderr,"Error parsing MD\n"); exit(0);
      }
      GT_PARSE_NUMBER(text_line,number);
      // Check CIGAR mismatches
      while (misms!=NULL && misms->position <= key_pos+number) { // Next
        // Account mismatch
        switch (misms->misms_type) {
          case MISMS: break;
          case INS:
            text_pos += misms->size;
            break;
          case DEL:
            key_pos += misms->size;
            break;
          default:
            break;
        }
        // Insert & next
        gt_vector_insert(misms_vector,*misms,gt_misms);
        gt_isp_parse_sam_opt_md_next_mismatch(map,&misms,&misms_pos);
      }
      // Next
      key_pos += number;
      text_pos += number;
    }
  }
  // Swap vectors & free
  GT_SWAP(map->mismatches,misms_vector);
  gt_vector_delete(misms_vector);
  return 0;
}
GT_INLINE gt_status gt_isp_parse_sam_opt_xa_bwa(
    const char** const text_line,gt_alignment* const alignment,
    gt_vector* const maps_vector,gt_sam_pending_end* const pending) {
  *text_line+=5;
  while (**text_line!=TAB && **text_line!=EOL) { // Read new attached maps
    gt_map* map = gt_map_new();
    map->phred_score = 0;
    map->gt_score = 0;
    gt_map_set_base_length(map,gt_alignment_get_read_length(alignment));
    // Sequence-name/Chromosome
    const char* const seq_name = *text_line;
    GT_READ_UNTIL(text_line,**text_line==COMA);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    gt_map_set_seq_name(map,seq_name,*text_line-seq_name);
    GT_NEXT_CHAR(text_line);
    // Position
    if (**text_line==MINUS) {
      gt_map_set_strand(map,REVERSE);
      GT_NEXT_CHAR(text_line);
    } else if (**text_line==PLUS) {
      gt_map_set_strand(map,FORWARD);
      GT_NEXT_CHAR(text_line);
    } else {
      gt_map_delete(map);
      return GT_ISP_PE_BAD_CHARACTER;
    }
    GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(map->position);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    GT_NEXT_CHAR(text_line);
    // CIGAR // TODO: Parse it !!
    GT_READ_UNTIL(text_line,**text_line==COMA);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    GT_NEXT_CHAR(text_line);
    // Edit distance
    GT_PARSE_NUMBER(text_line,map->gt_score);
    GT_READ_UNTIL(text_line,**text_line==SEMICOLON);
    if (**text_line==SEMICOLON) GT_NEXT_CHAR(text_line);
    // Add it to the list
    gt_vector_insert(maps_vector,map,gt_map*);
    // Consider pending relation
    ++pending->num_maps;
  }
  return 0;
}

GT_INLINE gt_status gt_isp_parse_sam_optional_field(
    const char** const text_line,gt_alignment* const alignment,
    gt_vector* const maps_vector,gt_sam_pending_end* const pending,
    const bool is_mapped) {
  const char* const init_opt_field = *text_line;

  /*
   * XA:Z:chr17,-34553512,125M,0;chr17,-34655077,125M,0;
   */
  GT_ISP_IF_OPT_FIELD(text_line,'X','A','Z') {
    if (!is_mapped) return GT_ISP_PE_SAM_UNMAPPED_XA;
    if (gt_isp_parse_sam_opt_xa_bwa(text_line,alignment,maps_vector,pending)) {
      *text_line = init_opt_field;
    }
  } GT_ISP_END_OPT_FIELD;

  GT_ISP_IF_OPT_FIELD(text_line,'M','D','Z') {
    if (is_mapped) {
      if (gt_isp_parse_sam_opt_md(text_line,alignment,maps_vector,pending)) {
        *text_line = init_opt_field;
      }
    }
  } GT_ISP_END_OPT_FIELD;

  /*
   * SAM-like field (store it as attribute and skip)
   */
  if (GT_INPUT_PARSER_IS_SAM_ATTRIBUTE(text_line)) {
    // Keep the beginning of the TAG
    const char* const attribute_start = *text_line;
    // Select proper attributes
    gt_attributes* attributes = NULL;
    // the attributes might not be initialized
    if (!is_mapped) {
        if(alignment->attributes==NULL){
          alignment->attributes = gt_attributes_new();
        }
        attributes = alignment->attributes;
    } else {
      gt_map* map = *gt_vector_get_last_elm(maps_vector,gt_map*);
      if(map->attributes==NULL){
        map->attributes = gt_attributes_new();
      }
      attributes = map->attributes;
    }
    // Parse OPT field
    if (!gt_input_parse_sam_optional_field(text_line,attributes)) return 0;
    *text_line = attribute_start;
  }

  // Skip the content as we cannot understand it (:-S)
  GT_READ_UNTIL(text_line,**text_line==TAB);

  return 0;
}

// TODO: Increase the level of checking SAM consistency
GT_INLINE gt_status gt_isp_parse_sam_alignment(
    char** const text_line,gt_template* const _template,gt_alignment* const _alignment,
    uint64_t* const alignment_flag,gt_sam_pending_end* const pending,const bool override_pairing) {
  gt_status error_code;
  bool is_mapped = true, is_single_segment;
  gt_map* map = gt_map_new();
  /*
   * Parse FLAG
   */
  if (!gt_is_number(**text_line)) {
    gt_map_delete(map);
    return GT_ISP_PE_EXPECTED_NUMBER;
  }
  GT_PARSE_NUMBER(text_line,*alignment_flag);
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  // Process flags
  const bool reverse_strand = (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT);
  is_mapped = !(*alignment_flag&GT_SAM_FLAG_UNMAPPED);
  is_single_segment = override_pairing || // TODO NovoAlign Doesn't need extra 3rd cnd
      !(*alignment_flag&GT_SAM_FLAG_MULTIPLE_SEGMENTS) || !(*alignment_flag&GT_SAM_FLAG_PROPERLY_ALIGNED);
  pending->end_position = (is_single_segment) ? 0 : ((*alignment_flag&GT_SAM_FLAG_FIRST_SEGMENT)?0:1);
  if (reverse_strand)  {
    gt_map_set_strand(map,REVERSE);
  } else {
    gt_map_set_strand(map,FORWARD);
  }
  // Allocate template/alignment handlers
  gt_alignment* alignment;
  if (_template) {
    alignment = gt_template_get_block_dyn(_template,0);
    if (pending->end_position==1) {
      alignment = gt_template_get_block_dyn(_template,1);
    }
  } else {
    GT_NULL_CHECK(_alignment);
    alignment = _alignment;
  }
  if (!gt_attributes_get(alignment->attributes,GT_ATTR_ID_SAM_FLAGS)) {
    gt_attributes_add(alignment->attributes,GT_ATTR_ID_SAM_FLAGS,alignment_flag,uint64_t);
  }
  /*
   * Parse RNAME (Sequence-name/Chromosome)
   */
  const char* const seq_name = *text_line;
  uint64_t seq_length = 0;
  if (gt_expect_false(**text_line==STAR)) {
    is_mapped=false; /* Unmapped */
    GT_NEXT_CHAR(text_line);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  } else {
    GT_READ_UNTIL(text_line,**text_line==TAB);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    seq_length = *text_line-seq_name;
    gt_map_set_seq_name(map,seq_name,seq_length);
    GT_NEXT_CHAR(text_line);
  }
  /*
   * Parse POS (Position 1-based)
   */
  GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(map->position);
  if (map->position==0) is_mapped=false; /* Unmapped */
  GT_NEXT_CHAR(text_line);
  /*
   * Parse MAPQ (Score)
   */
  GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(map->phred_score);
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  /*
   * Parse CIGAR
   */
  if ((error_code=gt_isp_parse_sam_cigar(text_line,&map,reverse_strand))) {
    gt_map_delete(map); return error_code;
  }
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  /*
   * Parse RNEXT (Sequence-name of the next segment)
   */
  if (**text_line==STAR || is_single_segment || !is_mapped ||
      (*alignment_flag&GT_SAM_FLAG_NEXT_UNMAPPED)) {
    gt_string_clear(&pending->next_seq_name);
    GT_ISP_PARSE_SAM_ALG_SKIP_FIELD(); // RNEXT
    GT_ISP_PARSE_SAM_ALG_SKIP_FIELD(); // PNEXT
  } else {
    // Parse RNEXT
    if (**text_line==EQUAL) {
      gt_string_set_nstring(&pending->next_seq_name,(char* const)seq_name,seq_length);
      GT_NEXT_CHAR(text_line);
      if (**text_line!=TAB) {
        gt_map_delete(map); return GT_ISP_PE_BAD_CHARACTER;
      }
      GT_NEXT_CHAR(text_line);
    } else {
      const char* const next_seq_name = *text_line;
      GT_READ_UNTIL(text_line,**text_line==TAB);
      GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
      gt_string_set_nstring(&pending->next_seq_name,(char* const)next_seq_name,*text_line-next_seq_name);
      GT_NEXT_CHAR(text_line);
    }
    // Parse PNEXT (Position of the next segment)
    GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(pending->next_position);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
    if (pending->next_position==0) {
      gt_string_clear(&pending->next_seq_name);
    } else {
      gt_string_set_nstring(&pending->map_seq_name,(char* const)seq_name,seq_length);
      pending->num_maps = 1;
      pending->map_position = gt_map_get_global_coordinate(map);
    }
  }
  /*
   * Parse TLEN (Template Length)
   */
  GT_ISP_PARSE_SAM_ALG_SKIP_FIELD();
  /*
   * Parse SEQ (READ)
   */
  if (gt_expect_false(**text_line==STAR)) {
    GT_NEXT_CHAR(text_line);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  } else {
    const char* const seq_read = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    const uint64_t read_length = *text_line-seq_read;
    GT_NEXT_CHAR(text_line);
    // Set the read
    if (gt_string_is_null(alignment->read)) {
      gt_string_set_nstring_static(alignment->read,seq_read,read_length);
      if (reverse_strand) {
        gt_dna_string_reverse_complement(alignment->read);
      }
    }
//    else if (pedantic) {
//      // Here we should check that the read if the same as previous occurrences
//      // Was disabled due to wrong SAM outputs produced some mappers when trimming
//      gt_string* const check_read = gt_string_new(0); // Writing to buffer's info
//      gt_string_set_nstring(check_read,seq_read,read_length);
//      if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT) {
//        gt_dna_string_reverse_complement(check_read);
//      }
//      const bool equals = gt_string_equals(alignment->read,check_read);
//      gt_string_delete(check_read);
//
//      if (!equals) {
//        gt_map_delete(map);
//        return GT_ISP_PE_WRONG_READ_CONTENT;
//      }
//    }
    if (gt_map_get_base_length(map)==0) gt_map_set_base_length(map,gt_alignment_get_read_length(alignment));
  }
  /*
   * Parse QUAL (QUALITY STRING)
   */
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
  if (gt_expect_false(**text_line==STAR && (*((*text_line)+1)==TAB || *((*text_line)+1)==EOL || *((*text_line)+1)==EOS) )) {
    GT_NEXT_CHAR(text_line);
  } else {
    const char* const seq_qual = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB);
    const uint64_t read_length = *text_line-seq_qual;
    if (gt_string_is_null(alignment->qualities)) {
      gt_string_set_nstring_static(alignment->qualities,seq_qual,read_length);
      if (reverse_strand) {
        gt_string_reverse(alignment->qualities);
      }
    }
    if (gt_map_get_base_length(map)==0) gt_map_set_base_length(map,gt_string_get_length(alignment->qualities));
  }
  if (!gt_string_is_null(alignment->read) && !gt_string_is_null(alignment->qualities)) {
    gt_fatal_check(gt_string_get_length(alignment->read)!=gt_string_get_length(alignment->qualities),ALIGNMENT_READ_QUAL_LENGTH);
  }
  // Build a list of alignments
  gt_vector *maps_vector = NULL;
  if (is_mapped) {
    maps_vector = gt_vector_new(10,sizeof(gt_map*));
    gt_vector_insert(maps_vector,map,gt_map*);
  } else {
    gt_map_delete(map);
  }
  /*
   * OPTIONAL FIELDS
   */
  while (**text_line==TAB) {
    GT_NEXT_CHAR(text_line);
    if ((error_code=gt_isp_parse_sam_optional_field(
        (const char** const)text_line,alignment,maps_vector,pending,is_mapped))) {
      if (maps_vector) {
        GT_VECTOR_ITERATE(maps_vector,map_elm,map_pos,gt_map*) gt_map_delete(*map_elm);
      }
      return error_code;
    }
  }
  // Add the main map
  if (**text_line!=EOL && **text_line!=EOS) return GT_ISP_PE_BAD_CHARACTER;
  if (is_mapped) {
    pending->map_displacement = gt_alignment_get_num_maps(alignment);
    if (override_pairing) {
      gt_alignment_insert_map_gt_vector(alignment,maps_vector);
    } else {
      GT_VECTOR_ITERATE(maps_vector,map_elm,map_pos,gt_map*) {
        gt_alignment_inc_counter(alignment,gt_map_get_global_distance(*map_elm));
        gt_alignment_add_map(alignment,*map_elm);
      }
    }
    gt_vector_delete(maps_vector);
  }
  /*
   * Reverse CIGAR (if needed)
   *   Just before XA due to MD synch
   */
  if (reverse_strand) {
    GT_MAP_ITERATE(map,map_it) {
      gt_map_reverse_misms(map_it);
    }
  }
  return 0;
}

GT_INLINE bool gt_isp_fetch_next_line(
    gt_buffered_input_file* const buffered_sam_input,gt_string* const expected_tag,const bool chomp_tag) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_NULL_CHECK(expected_tag);
  // Check next record/line
  gt_input_sam_parser_next_record(buffered_sam_input);
  if (gt_buffered_input_file_eob(buffered_sam_input)) return false;
  // Fetch next tag
  gt_string* const next_tag = gt_string_new(16); // initialize so it not static ?
  char* ptext_line;
  if (gt_isp_read_tag((const char** const)&(buffered_sam_input->cursor),
      (const char** const)&ptext_line,next_tag)) return false;
  if (chomp_tag) gt_input_parse_tag_chomp_pairend_info(next_tag);
  const bool same_tag = gt_string_equals(expected_tag,next_tag);
  gt_string_delete(next_tag);
  if (same_tag) {
    buffered_sam_input->cursor = ptext_line;
    return true;
  } else {
    return false;
  }
}

GT_INLINE void gt_isp_add_mmap(
    gt_template* const template,const uint64_t start_pos_end1,const uint64_t start_pos_end2,
    const uint64_t pending_maps_end1,const uint64_t pending_maps_end2) {
  gt_map** mmap_end1 = gt_vector_get_elm(gt_template_get_block_dyn(template,0)->maps,start_pos_end1,gt_map*);
  gt_map** mmap_end2 = gt_vector_get_elm(gt_template_get_block_dyn(template,1)->maps,start_pos_end2,gt_map*);
  uint64_t i;
  gt_mmap_attributes attr;
  const uint64_t pending_maps = GT_MAX(pending_maps_end1,pending_maps_end2);
  for (i=0;i<pending_maps;++i) {
    gt_map* map_end[2];
    map_end[0] = (pending_maps_end1>1) ? mmap_end1[i] : mmap_end1[0];
    map_end[1] = (pending_maps_end2>1) ? mmap_end2[i] : mmap_end2[0];
    attr.distance = gt_map_get_global_distance(map_end[0])+gt_map_get_global_distance(map_end[1]);
    attr.phred_score = map_end[0]->phred_score;
    gt_template_inc_counter(template,attr.distance);
    gt_template_add_mmap_ends(template,map_end[0],map_end[1],&attr);
  }
}

GT_INLINE bool gt_isp_check_pending_record__add_mmap(
    gt_template* const template,gt_sam_pending_end* const pending,
    const uint64_t end_position,gt_string* const seq_name,const uint64_t position,
    const uint64_t map_displacement,const uint64_t num_maps) {
  if (pending->end_position!=end_position &&
      gt_string_equals(&pending->next_seq_name,seq_name) &&
      pending->next_position==position) { // Found!
    // BWA_Compact. MMAPs paired against MMaps need to have the same cardinality (otherwise it's unpaired)
    if (pending->num_maps!=1 && num_maps!=1 && pending->num_maps!=num_maps) return false;
    // Insert mmap(s)
    if (pending->end_position==1) {
      gt_isp_add_mmap(template,map_displacement,pending->map_displacement,num_maps,pending->num_maps);
    } else {
      gt_isp_add_mmap(template,pending->map_displacement,map_displacement,pending->num_maps,num_maps);
    }
    return true;
  }
  return false;
}

GT_INLINE void gt_isp_solve_pending_maps(
    gt_vector* pending_v,gt_sam_pending_end* pending,gt_template* const template) {
  bool found_match = false;
  // Look into pending records
  GT_VECTOR_ITERATE(pending_v,pending_elm,pending_counter,gt_sam_pending_end) {
    if (gt_string_is_null(&pending_elm->next_seq_name)) continue;
    if ((found_match=gt_isp_check_pending_record__add_mmap(
            template,pending_elm,pending->end_position,&pending->map_seq_name,
            pending->map_position,pending->map_displacement,pending->num_maps))) {
      gt_string_clear(&pending_elm->next_seq_name); // Mark as solved
      break;
    }
  }
  // Queue if not found
  if (!found_match) gt_vector_insert(pending_v,*pending,gt_sam_pending_end);
}

GT_INLINE gt_status gt_isp_solve_remaining_maps(gt_vector* const pending_v,gt_template* const template) {
  gt_status error_code = 0;
  GT_VECTOR_ITERATE(pending_v,pending_elm,pending_counter,gt_sam_pending_end) {
    if (!gt_string_is_null(&pending_elm->next_seq_name)) {
      // Look the pending map in the already stored maps (BWA-Based)
      const uint64_t map_end = (pending_elm->end_position+1)%2;
      uint64_t pos = 0;
      bool found = false;
      GT_ALIGNMENT_ITERATE(gt_template_get_block_dyn(template,map_end),map) {
        if ((found=gt_isp_check_pending_record__add_mmap(template,pending_elm,map_end,
            gt_map_get_string_seq_name(map),gt_map_get_global_coordinate(map),pos,1))) break;
        ++pos;
      }
      // Check solved
      if (!found) {
        error_code = GT_ISP_PE_UNSOLVED_PENDING_MAPS;
        break;
      }
    }
  }
  return error_code;
}

#define gt_isp_skip_remaining_records(buffered_sam_input,tag) while (gt_isp_fetch_next_line(buffered_sam_input,tag,false))

/* SAM general */
GT_INLINE gt_status gt_input_sam_parser_parse_template(
    gt_buffered_input_file* const buffered_sam_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  gt_status error_code;
  char** const text_line = (char** const)&(buffered_sam_input->cursor);
  // Read initial TAG (QNAME := Query template)
  if ((error_code=gt_isp_read_tag((const char** const)text_line,(const char** const)text_line,template->tag))) return error_code;
  gt_input_parse_tag_chomp_pairend_info(template->tag);
  // Read all maps related to this TAG
  gt_vector* pending_v = gt_vector_new(GT_ISP_NUM_INITIAL_MAPS,sizeof(gt_sam_pending_end));
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending = GT_SAM_INIT_PENDING;
    uint64_t alignment_flag;
    if (gt_expect_false(error_code=gt_isp_parse_sam_alignment(
          text_line,template,NULL,&alignment_flag,&pending,false))) {
      gt_vector_delete(pending_v);
      gt_isp_skip_remaining_records(buffered_sam_input,template->tag);
      return error_code;
    }
    // Solve pending ends
    if (!gt_string_is_null(&pending.next_seq_name)) gt_isp_solve_pending_maps(pending_v,&pending,template);
  } while (gt_isp_fetch_next_line(buffered_sam_input,template->tag,true));
  // Check for unsolved pending maps (try to solve them)
  error_code = gt_isp_solve_remaining_maps(pending_v,template);
  gt_vector_delete(pending_v);
  if (error_code) gt_isp_skip_remaining_records(buffered_sam_input,template->tag);
  // Setup alignment's tag info
  gt_template_setup_pair_attributes_to_alignments(template,true);
  return error_code;
}
/* SOAP2-SAM */
GT_INLINE gt_status gt_input_sam_parser_parse_soap_template(
    gt_buffered_input_file* const buffered_sam_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  char** const text_line = (char** const)&(buffered_sam_input->cursor);
  gt_status error_code;
  // Read initial TAG (QNAME := Query template)
  if ((error_code=gt_isp_read_tag((const char** const)text_line,(const char** const)text_line,template->tag))) return error_code;
  gt_input_parse_tag_chomp_pairend_info(template->tag);
  // Read all maps related to this TAG
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending = GT_SAM_INIT_PENDING;
    uint64_t alignment_flag;
    if (gt_expect_false(error_code=gt_isp_parse_sam_alignment(
          text_line,template,NULL,&alignment_flag,&pending,false))) {
      gt_isp_skip_remaining_records(buffered_sam_input,template->tag);
      return error_code;
    }
  } while (gt_isp_fetch_next_line(buffered_sam_input,template->tag,true));
  // SOAP2 paired maps convention. Add maps
  gt_alignment* const alignment_end0 = gt_template_get_block_dyn(template,0);
  gt_alignment* const alignment_end1 = gt_template_get_block_dyn(template,1);
  if (gt_alignment_get_num_maps(alignment_end0) !=
      gt_alignment_get_num_maps(alignment_end1)) return GT_ISP_PE_UNSOLVED_PENDING_MAPS;
  uint64_t pos_end_it=0;
  GT_ALIGNMENT_ITERATE(alignment_end0,map_end1) {
    gt_map* map_end2 = gt_alignment_get_map(alignment_end1,pos_end_it);
    gt_mmap_attributes attr;
    attr.distance = gt_map_get_global_distance(map_end1)+gt_map_get_global_distance(map_end2);
    gt_template_add_mmap_ends(template,map_end1,map_end2,&attr);
    ++pos_end_it;
  }
  // Setup alignment's tag info
  gt_template_setup_pair_attributes_to_alignments(template,true);
  return 0;
}
/* SE-SAM */
GT_INLINE gt_status gt_input_sam_parser_parse_alignment(
    gt_buffered_input_file* const buffered_sam_input,gt_alignment* alignment) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_ALIGNMENT_CHECK(alignment);
  gt_status error_code;
  char** const text_line = (char** const)&(buffered_sam_input->cursor);
  // Read initial TAG (QNAME := Query template)
  if ((error_code=gt_isp_read_tag((const char** const)text_line,(const char** const)text_line,alignment->tag))) return error_code;
  // Read all maps related to this TAG
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending = GT_SAM_INIT_PENDING;
    uint64_t alignment_flag;
    if (gt_expect_false((error_code=gt_isp_parse_sam_alignment(
        text_line,NULL,alignment,&alignment_flag,&pending,true))!=0)) {
      return error_code;
    }
  } while (gt_isp_fetch_next_line(buffered_sam_input,alignment->tag,false));
  // Chomp /1/2 and add the pair info
  int64_t pair = gt_input_parse_tag_chomp_pairend_info(alignment->tag);
  if (pair) gt_attributes_add(alignment->attributes,GT_ATTR_ID_TAG_PAIR,&pair,int64_t);
  return 0;
}


/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_sam_parser_get_template(
    gt_buffered_input_file* const buffered_sam_input,gt_template* const template,gt_sam_parser_attributes* const sam_parser_attr) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(sam_parser_attr);
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_sam_input)) {
    if ((error_code=gt_input_sam_parser_reload_buffer(buffered_sam_input))!=GT_ISP_OK) return error_code;
  }
  // Check file format
  gt_input_file* input_file = buffered_sam_input->input_file;
  if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
    gt_error(PARSE_SAM_BAD_FILE_FORMAT,gt_input_file_get_file_name(input_file),buffered_sam_input->current_line_num,(uint64_t)0);
    return GT_ISP_FAIL;
  }
  // Prepare the template
  char* const line_start = buffered_sam_input->cursor;
  const uint64_t line_num = buffered_sam_input->current_line_num;
  gt_template_clear(template,true);
  template->template_id = line_num;
  // Parse template
  if (gt_expect_false(sam_parser_attr->sam_soap_style)) {
    error_code=gt_input_sam_parser_parse_soap_template(buffered_sam_input,template);
  } else {
    error_code=gt_input_sam_parser_parse_template(buffered_sam_input,template);
  }
  if (error_code) {
    gt_input_sam_parser_prompt_error(buffered_sam_input,line_num,
        buffered_sam_input->cursor-line_start,error_code);
    gt_input_sam_parser_next_record(buffered_sam_input);
    return GT_ISP_FAIL;
  }
  return GT_ISP_OK;
}

GT_INLINE gt_status gt_input_sam_parser_get_alignment(
    gt_buffered_input_file* const buffered_sam_input,gt_alignment* const alignment,gt_sam_parser_attributes* const sam_parser_attr) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_ALIGNMENT_CHECK(alignment);
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_sam_input)) {
    if ((error_code=gt_input_sam_parser_reload_buffer(buffered_sam_input))!=GT_ISP_OK) return error_code;
  }
  // Check file format
  gt_input_file* input_file = buffered_sam_input->input_file;
  if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
    gt_error(PARSE_SAM_BAD_FILE_FORMAT,gt_input_file_get_file_name(input_file),buffered_sam_input->current_line_num,(uint64_t)0);
    return GT_ISP_FAIL;
  }
  // Allocate memory for the alignment
  char* const line_start = buffered_sam_input->cursor;
  const uint64_t line_num = buffered_sam_input->current_line_num;
  gt_alignment_clear(alignment);
  alignment->alignment_id = line_num;
  // Parse alignment
  if ((error_code=gt_input_sam_parser_parse_alignment(buffered_sam_input,alignment))) {
    gt_input_sam_parser_prompt_error(buffered_sam_input,line_num,
        buffered_sam_input->cursor-line_start,error_code);
    gt_input_sam_parser_next_record(buffered_sam_input);
    return GT_ISP_FAIL;
  }
  return GT_ISP_OK;
}

