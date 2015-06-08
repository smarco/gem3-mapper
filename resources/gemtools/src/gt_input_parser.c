/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_tag_parser.c
 * DATE: 01/06/2012
 * AUTHOR(S): Thasso Griebel <thasso.griebel@gmail.com>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_input_parser.h"
#include "gt_input_fasta_parser.h"

/*
 * Basic Parsing Functions
 */
GT_INLINE void gt_input_parse_next_char(const char** const text_line) {
  GT_NULL_CHECK(text_line);
  GT_NEXT_CHAR(text_line);
}
GT_INLINE void gt_input_parse_skip_chars(const char** const text_line,uint64_t num_chars) {
  GT_NULL_CHECK(text_line);
  while ((num_chars--) > 0 && !GT_IS_EOL(text_line)) GT_NEXT_CHAR(text_line);
}
GT_INLINE gt_status gt_input_parse_eol(const char** const text_line) {
  GT_NULL_CHECK(text_line);
  return GT_IS_EOL(text_line);
}
GT_INLINE void gt_input_parse_field(const char** const text_line,const char delimiter,gt_string* const string) {
  GT_NULL_CHECK(text_line);
  // Read field
  const char* const string_begin = *text_line;
  while (gt_expect_true(**text_line!=delimiter && !GT_IS_EOL(text_line))) GT_NEXT_CHAR(text_line);
  // Copy string
  if (string) gt_string_set_nstring_static(string,string_begin,(*text_line-string_begin));
  // Skip delimiter
  if (**text_line==delimiter) GT_NEXT_CHAR(text_line);
}
GT_INLINE gt_status gt_input_parse_integer(const char** const text_line,int64_t* const value) {
  GT_NULL_CHECK(text_line);
  GT_NULL_CHECK(value);
  int64_t number = 0;
  if (**text_line=='0' && (*(*text_line+1)=='x' || *(*text_line+1)=='X')) {
    *text_line+=2;
    if (gt_expect_false(!gt_is_hex_digit(**text_line))) return -1;
    // Parse Value
    while (gt_expect_true(gt_is_hex_digit(**text_line))) {
      number = (number<<4) + gt_get_hex_cipher(**text_line);
      GT_NEXT_CHAR(text_line);
    }
  } else {
    // Parse sign (if any)
    bool positive = true;
    if (gt_expect_false(**text_line==PLUS)) {
      GT_NEXT_CHAR(text_line);
    } else if (gt_expect_false(**text_line==MINUS)) {
      positive = false;
      GT_NEXT_CHAR(text_line);
    }
    // Parse number
    if (gt_expect_false(!gt_is_number(**text_line))) return -1;
    while (gt_expect_true(gt_is_number(**text_line))) {
      number = (number*10) + gt_get_cipher(**text_line);
      GT_NEXT_CHAR(text_line);
    }
    // Add sign
    if (gt_expect_false(!positive)) number = -number;
  }
  *value = number;
  return 0;
}
GT_INLINE gt_status gt_input_parse_double(const char** const text_line,double* const value) {
  GT_NULL_CHECK(text_line);
  GT_NULL_CHECK(value);
  /*
   * [+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)
   *   Sign ::= [+-]
   *   Integer ::= [0-9]+
   *   Dot ::= "."
   *   Decimal ::= [0-9]+
   *   Exponent ::= [0-9]+
   */
  // Parse sign (if any)
  bool positive = true;
  if (gt_expect_false(**text_line==PLUS)) {
    GT_NEXT_CHAR(text_line);
  } else if (gt_expect_false(**text_line==MINUS)) {
    positive = false;
    GT_NEXT_CHAR(text_line);
  }
  // Parse Integer
  double integer = 0.0;
  if (gt_expect_true(gt_is_number(**text_line))) { // 45
    while (gt_expect_true(gt_is_number(**text_line))) {
      integer = (integer*10) + gt_get_cipher(**text_line);
      GT_NEXT_CHAR(text_line);
    }
  }
  // Dot
  if (gt_expect_true(**text_line==DOT)) { // .045
    GT_NEXT_CHAR(text_line);
    // Dot forces to have a decimal part
    if (gt_expect_false(!gt_is_number(**text_line))) return -1;
  }
  // Parse decimal
  double decimal = 0.0;
  if (gt_expect_true(gt_is_number(**text_line))) { // 45
    while (gt_expect_true(gt_is_number(**text_line))) {
      decimal = (decimal/10) + gt_get_cipher(**text_line);
      GT_NEXT_CHAR(text_line);
    }
  }
  // Parse exponent
  double exponent = 0.0;
  if (gt_expect_false(**text_line=='e' || **text_line=='E')) {
    GT_NEXT_CHAR(text_line);
    // Parse sign (if any)
    bool exponent_positive = true;
    if (gt_expect_false(**text_line==PLUS)) {
      GT_NEXT_CHAR(text_line);
    } else if (gt_expect_false(**text_line==MINUS)) {
      exponent_positive = false;
      GT_NEXT_CHAR(text_line);
    }
    // Parse number
    if (gt_expect_false(!gt_is_number(**text_line))) return -1;
    while (gt_expect_true(gt_is_number(**text_line))) {
      exponent = (exponent*10) + gt_get_cipher(**text_line);
      GT_NEXT_CHAR(text_line);
    }
    // Apply sign
    if (!exponent_positive) exponent = -exponent;
  }
  // Compose number
  *value = integer + decimal;
  if (exponent!=0) *value = *value * pow(10,exponent);
  if (!positive) *value = -(*value);
  return 0;
}
/*
 * Count number of ':' in field + 1
 */
GT_INLINE uint64_t gt_input_parse_count_colons_in_field(const char* const text_line) {
  uint64_t i = 0, count = 0;
  while (text_line[i]!=TAB && text_line[i]!=SPACE && text_line[i]!=EOL && text_line[i]!=EOS) {
    if (text_line[i]==':') ++count;
    ++i;
  }
  return count;
}
/*
 * Read trim. Covers
 *   lt:Z:5:ACGTA:BBBBB lt:Z:5:ACGTA: lt:Z:5::BBBBB lt:Z:5:: lt:Z:5
 *   rt:Z:5:ACGTA:BBBBB rt:Z:5:ACGTA: rt:Z:5::BBBBB rt:Z:5:: rt:Z:5
 */
GT_INLINE gt_status gt_input_parse_attribute_trim(const char** const text_line,gt_read_trim* const trim_info) {
  GT_NULL_CHECK(text_line);
  GT_NULL_CHECK(trim_info);
  (*text_line)+=5; // Skip lt:Z: | rt:Z:
  /*
   * Parse length
   */
  GT_PARSE_NUMBER(text_line,trim_info->length);
  if (GT_IS_EOL(text_line) || **text_line==SPACE || **text_line==TAB) return 0;
  /*
   * Parse trimmed read
   */
  if (**text_line!=COLON) return GT_PE_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  const char* const trimmed_read_begin = *text_line;
  GT_READ_UNTIL(text_line,**text_line==COLON || **text_line==SPACE || **text_line==TAB);
  const uint64_t trimmed_read_length = *text_line-trimmed_read_begin;
  if (trimmed_read_length!=trim_info->length) return GT_PE_BAD_TRIM_READ_STRING_LENGTH;
  // Set trimmed read
  trim_info->trimmed_read = gt_string_new(trimmed_read_length+1);
  gt_string_set_nstring_static(trim_info->trimmed_read,trimmed_read_begin,trimmed_read_length);
  if (GT_IS_EOL(text_line) || **text_line==SPACE || **text_line==TAB) return 0;
  /*
   * Parse trimmed qualities
   */
  if (**text_line!=COLON) return GT_PE_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  const char* const trimmed_qual_begin = *text_line;
  GT_READ_UNTIL(text_line,**text_line==SPACE || **text_line==TAB);
  const uint64_t trimmed_qual_length = *text_line-trimmed_qual_begin;
  if (trimmed_qual_length!=trim_info->length) return GT_PE_BAD_TRIM_QUAL_STRING_LENGTH;
  // Set trimmed qualities
  trim_info->trimmed_qualities = gt_string_new(trimmed_qual_length+1);
  gt_string_set_nstring_static(trim_info->trimmed_qualities,trimmed_qual_begin,trimmed_qual_length);
  if (GT_IS_EOL(text_line) || **text_line==SPACE || **text_line==TAB) {
    return 0;
  } else {
    return GT_PE_BAD_CHARACTER;
  }
}
/*
 * Read Chunks. Covers
 *   sr:Z:1:3
 */
GT_INLINE gt_status gt_input_parse_attribute_segmented_read(const char** const text_line,gt_segmented_read_info* const segmented_read_info) {
  GT_NULL_CHECK(text_line);
  GT_NULL_CHECK(segmented_read_info);
  (*text_line)+=5; // Skip sr:Z:
  /*
   * Parse Segment ID
   */
  GT_PARSE_NUMBER(text_line,segmented_read_info->segment_id);
  if (GT_IS_EOL(text_line) || **text_line==SPACE || **text_line==TAB) return 0;
  /*
   * Parse total number of Segments
   */
  if (**text_line!=COLON) return GT_PE_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  GT_PARSE_NUMBER(text_line,segmented_read_info->total_segments);
  if (GT_IS_EOL(text_line) || **text_line==SPACE || **text_line==TAB) {
    return 0;
  } else {
    return GT_PE_BAD_CHARACTER;
  }
}
/*
 * Parse any SAM-like attribute
 */
GT_INLINE gt_status gt_input_parse_sam_optional_field(const char** const text_line,gt_attributes* const attributes) {
  GT_NULL_CHECK(text_line);
  GT_ATTRIBUTES_CHECK(attributes);
  // Parse TAG
  const char* const tag = *text_line;
  (*text_line)+=2;
  // COLON
  if (**text_line!=COLON) return GT_PE_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  // Parse ID
  const char type_id = **text_line;
  GT_NEXT_CHAR(text_line);
  // COLON
  if (**text_line!=COLON) return GT_PE_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  // Parse Value
  switch (type_id) {
    case 'i': { // Integer value
      int64_t value;
      if (gt_input_parse_integer(text_line,&value)) return GT_PE_BAD_CHARACTER;
      gt_attributes_add_sam_ivalue(attributes,tag,type_id,value);
      break;
    }
    case 'f': { // Float value
      double value;
      if (gt_input_parse_double(text_line,&value)) return GT_PE_BAD_CHARACTER;
      gt_attributes_add_sam_fvalue(attributes,tag,type_id,value);
      break;
    }
    default: { // Otherwise is a string
      const char* const begin_value = *text_line;
      GT_READ_UNTIL(text_line,**text_line==SPACE || **text_line==TAB);
      const uint64_t length = (*text_line-begin_value);
      gt_string* const string = gt_string_new(length+1);
      gt_string_set_nstring_static(string,begin_value,length);
      gt_attributes_add_sam_svalue(attributes,tag,type_id,string);
      break;
    }
  }
  return 0;
}
/*
 * GT tag parsing
 *   Removes all pair info /1/2/3 and puts it as attribute
 *   Parses all TAG attributes
 *     Attributes ...
 *     CASAVA information
 *     Extra information ...
 */
#define GT_INPUT_PARSE_SAM_ATTRIBUTE_FIELD(text_line,char1,char2,type_char) \
  if ((*text_line)[0]==char1 && (*text_line)[1]==char2 && (*text_line)[2]==COLON && (*text_line)[3]==type_char && (*text_line)[4]==COLON)
GT_INLINE gt_status gt_input_parse_tag(const char** const text_line,gt_string* const tag,gt_attributes* const attributes) {
  GT_NULL_CHECK(text_line);
  GT_STRING_CHECK(tag);
  GT_ATTRIBUTES_CHECK(attributes);
  // Delimit the tag
  register uint64_t i = 0;
  const char* const tag_begin = *text_line;
  // Parse Tag
  GT_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE); // Read until first SPACE or TAB
  const uint64_t tag_length = *text_line-tag_begin;
  gt_string_set_nstring_static(tag,tag_begin,tag_length);
  // Add pair info and chomp /1/2/3 info (if any)
  const int64_t tag_pair = gt_input_parse_tag_chomp_pairend_info(tag);
  gt_attributes_add(attributes,GT_ATTR_ID_TAG_PAIR,&tag_pair,int64_t);
  gt_string_append_eos(tag);
  /*
   * Parse all extra TAG-info
   */
  while (!GT_IS_EOL(text_line) && **text_line!=TAB) {
    GT_NEXT_CHAR(text_line); // Skip space
    // LEFT-Trim
    GT_INPUT_PARSE_SAM_ATTRIBUTE_FIELD(text_line,'l','t','Z') {
      const char* const attribute_start = *text_line;
      gt_read_trim trim_info;
      if (gt_input_parse_attribute_trim(text_line,&trim_info)) {
        *text_line = attribute_start;
      } else {
        if (trim_info.length>0) gt_attributes_annotate_left_trim(attributes,&trim_info);
        continue;
      }
    }
    // RIGHT-Trim
    GT_INPUT_PARSE_SAM_ATTRIBUTE_FIELD(text_line,'r','t','Z') {
      const char* const attribute_start = *text_line;
      gt_read_trim trim_info;
      if (gt_input_parse_attribute_trim(text_line,&trim_info)) {
        *text_line = attribute_start;
      } else {
        if (trim_info.length>0) gt_attributes_annotate_right_trim(attributes,&trim_info);
        continue;
      }
    }
    // Read Chunks
    GT_INPUT_PARSE_SAM_ATTRIBUTE_FIELD(text_line,'s','r','Z') {
      const char* const attribute_start = *text_line;
      gt_segmented_read_info segmented_read_info;
      if (gt_input_parse_attribute_segmented_read(text_line,&segmented_read_info)) {
        *text_line = attribute_start;
      } else {
        gt_attributes_add(attributes,GT_ATTR_ID_SEGMENTED_READ_INFO,&segmented_read_info,gt_segmented_read_info);
        continue;
      }
    }
    /*
     * Any SAM-like attribute
     */
    if (GT_INPUT_PARSER_IS_SAM_ATTRIBUTE(text_line)) {
      const char* const attribute_start = *text_line;
      if (!gt_input_parse_sam_optional_field(text_line,attributes)) continue;
      *text_line = attribute_start;
    }
    /*
     * CASAVA 1.8 Attributes. @EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACG
     */
    if (gt_input_parse_count_colons_in_field(*text_line)==3) {
      const char* const casava_info_begin = *text_line;
      const int64_t pair = (casava_info_begin[0]=='1') ? GT_PAIR_PE_1 :
        ((casava_info_begin[0]=='2' || casava_info_begin[0]=='3') ? GT_PAIR_PE_2 : GT_PAIR_SE);
      if (pair==GT_PAIR_PE_1 || pair==GT_PAIR_PE_2) {
        GT_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE);
        const uint64_t casava_info_length = *text_line-casava_info_begin;
        gt_string* const casava_string = gt_string_new(casava_info_length+1);
        gt_string_set_nstring_static(casava_string,casava_info_begin,casava_info_length);
        gt_attributes_add_string(attributes,GT_ATTR_ID_TAG_CASAVA,casava_string);
        gt_attributes_add(attributes,GT_ATTR_ID_TAG_PAIR,&pair,int64_t);
        continue; // Next!
      }
    }
    /*
     * TAG Extra
     *   GT-41 add additional check to see if any extra attributes end in /1 /2 /3
     *   and no pair found yet if thats the case, this takes over, /1/2/3 is cut away
     *   and the tag is reset to the original tag but spaces are replaced with _
     *   E.g
     *     @SRR384920.1 HWI-ST382_0049:1:1:1217:1879/1
     *     @SRR384920.1 HWI-ST382_0049:1:1:1217:1879/2
     */
    const char* const extra_tag_begin = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE);
    const uint64_t extra_tag_length = *text_line-extra_tag_begin;
    gt_string* const extra_string = gt_string_new(extra_tag_length+1);
    gt_string_set_nstring_static(extra_string,extra_tag_begin,extra_tag_length);
    // Process extra tag information
    const int64_t tag_extra_pair = gt_input_parse_tag_chomp_pairend_info(extra_string);
    if (tag_extra_pair==GT_PAIR_PE_1 || tag_extra_pair==GT_PAIR_PE_2) {
      // Append to the tag an '_' plus the extra information
      gt_string_append_char(tag,UNDERSCORE);
      gt_string_append_gt_string(tag,extra_string);
      // Replace all spaces
      const uint64_t tag_length = gt_string_get_length(tag);
      for (i=0;i<tag_length;i++) {
        if (tag->buffer[i]==SPACE) tag->buffer[i] = UNDERSCORE;
      }
      gt_string_append_eos(tag);
      // Set pair info
      gt_attributes_add(attributes,GT_ATTR_ID_TAG_PAIR,&tag_extra_pair,int64_t);
      // Free
      gt_string_delete(extra_string);
    } else {
      gt_string* const attribute_extra_string = gt_attributes_get(attributes,GT_ATTR_ID_TAG_EXTRA);
      if (attribute_extra_string==NULL) {
        gt_attributes_add_string(attributes,GT_ATTR_ID_TAG_EXTRA,extra_string);
      } else {
        gt_string_append_char(attribute_extra_string,SPACE);
        gt_string_append_gt_string(attribute_extra_string,extra_string);
        gt_string_delete(extra_string);
      }
    }
  } /* while (not end of tags) */
  GT_NEXT_CHAR(text_line);
  return GT_STATUS_OK; // Return OK
}
/*
 * Parse the end information {/1,/2}
 */
GT_INLINE uint64_t gt_input_parse_tag_chomp_pairend_info(gt_string* const tag) {
  GT_STRING_CHECK(tag);
  const uint64_t tag_length = gt_string_get_length(tag);
  if (tag_length>2 && *gt_string_char_at(tag,tag_length-2)==SLASH) {
    const char tag_end = *gt_string_char_at(tag,tag_length-1);
    if (tag_end=='1') {
      gt_string_set_length(tag,tag_length-2);
      return GT_PAIR_PE_1;
    } else if (tag_end=='2' || tag_end=='3') {
      gt_string_set_length(tag,tag_length-2);
      return GT_PAIR_PE_2;
    } else {
      return GT_PAIR_SE;
    }
  } else {
    return GT_PAIR_SE;
  }
}
