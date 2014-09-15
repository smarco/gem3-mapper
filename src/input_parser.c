/*
 * PROJECT: GEMMapper
 * FILE: input_parser.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "input_parser.h"

/*
 * Basic InputFile-Parsing Functions
 */
GEM_INLINE bool input_file_parse_next_char(input_file_t* const input_file) {
  INPUT_FILE_CHECK(input_file);
  return input_file_next_char(input_file);
}
GEM_INLINE error_code_t input_file_parse_skip_separators(input_file_t* const input_file) {
  INPUT_FILE_CHECK(input_file);
  if (gem_expect_false(input_file_eof(input_file))) return -1;
  const char current_char = input_file_get_current_char(input_file);
  if (gem_expect_false(current_char!=TAB && current_char!=SPACE)) return -1;
  // Skip until different from TAB or SPACE
  while (input_file_next_char(input_file)) {
    // Read character
    const char current_char = input_file_get_current_char(input_file);
    if (IS_ANY_EOL(current_char) || (current_char!=TAB && current_char!=SPACE)) break;
  }
  return 0;
}
GEM_INLINE void input_file_parse_skip_chars(input_file_t* const input_file,uint64_t num_chars) {
  INPUT_FILE_CHECK(input_file);
  while (!input_file_eof(input_file) && num_chars>0) {
    // Read character
    const char current_char = input_file_get_current_char(input_file);
    if (IS_ANY_EOL(current_char)) break;
    // Next
    input_file_next_char(input_file);
    --num_chars;
  }
}
GEM_INLINE void input_file_parse_skip_line(input_file_t* const input_file) {
  INPUT_FILE_CHECK(input_file);
  while (!input_file_eof(input_file)) {
    // Read character
    const char current_char = input_file_get_current_char(input_file);
    if (IS_ANY_EOL(current_char)) break;
    // Next
    input_file_next_char(input_file);
  }
}
GEM_INLINE bool input_file_parse_is_eol(input_file_t* const input_file) {
  INPUT_FILE_CHECK(input_file);
  return input_file_eof(input_file) || IS_ANY_EOL(input_file_get_current_char(input_file));
}
GEM_INLINE void input_file_parse_field(input_file_t* const input_file,const char delimiter,string_t* const string) {
  INPUT_FILE_CHECK(input_file);
  while (!input_file_eof(input_file)) {
    // Read character
    const char current_char = input_file_get_current_char(input_file);
    if (IS_ANY_EOL(current_char) || current_char==delimiter) break;
    // Append character
    string_append_char(string,current_char);
    // Next
    input_file_next_char(input_file);
  }
  string_append_eos(string);
}
GEM_INLINE error_code_t input_file_parse_integer(input_file_t* const input_file,int64_t* const value) {
  INPUT_FILE_CHECK(input_file);
  GEM_CHECK_NULL(value);
  int64_t number = 0;
  if (input_file_eof(input_file)) return -1;
  char current_char = input_file_get_current_char(input_file);
  // Parse sign (if any)
  bool positive = true;
  if (gem_expect_false(current_char==PLUS)) {
    if (!input_file_next_char(input_file)) return -1;
    current_char = input_file_get_current_char(input_file);
  } else if (gem_expect_false(current_char==MINUS)) {
    positive = false;
    if (!input_file_next_char(input_file)) return -1;
    current_char = input_file_get_current_char(input_file);
  }
  // Parse number
  if (gem_expect_false(!IS_DIGIT(current_char))) return -1;
  while (gem_expect_true(IS_DIGIT(current_char))) {
    number = (number*10) + GET_DIGIT(current_char);
    if (!input_file_next_char(input_file)) break;
    current_char = input_file_get_current_char(input_file);
  }
  // Add sign
  if (gem_expect_false(!positive)) number = -number;
  *value = number;
  return 0;
}

/*
 * Basic Text-Parsing Functions
 */
GEM_INLINE void input_text_parse_next_char(const char** const text_line) {
  GEM_CHECK_NULL(text_line);
  PARSER_NEXT_CHAR(text_line);
}
GEM_INLINE void input_text_parse_skip_chars(const char** const text_line,uint64_t num_chars) {
  GEM_CHECK_NULL(text_line);
  while ((num_chars--) > 0 && !PARSER_IS_EOL(text_line)) PARSER_NEXT_CHAR(text_line);
}
GEM_INLINE void input_text_parse_skip_line(const char** const text_line) {
  GEM_CHECK_NULL(text_line);
  PARSER_SKIP_LINE(text_line);
}
GEM_INLINE bool input_text_parse_is_eol(const char** const text_line) {
  GEM_CHECK_NULL(text_line);
  return PARSER_IS_EOL(text_line);
}
GEM_INLINE void input_text_parse_field(const char** const text_line,const char delimiter,string_t* const string) {
  GEM_CHECK_NULL(text_line);
  // Read field
  const char* const string_begin = *text_line;
  while (gem_expect_true(**text_line!=delimiter && !PARSER_IS_EOL(text_line))) PARSER_NEXT_CHAR(text_line);
  // Copy string
  if (string) string_set_buffer_const(string,string_begin,(*text_line-string_begin));
  // Skip delimiter
  if (**text_line==delimiter) PARSER_NEXT_CHAR(text_line);
}
GEM_INLINE error_code_t input_text_parse_integer(const char** const text_line,int64_t* const value) {
  GEM_CHECK_NULL(text_line);
  GEM_CHECK_NULL(value);
  int64_t number = 0;
  if (**text_line=='0' && (*(*text_line+1)=='x' || *(*text_line+1)=='X')) {
    *text_line+=2;
    if (gem_expect_false(!IS_HEX_DIGIT(**text_line))) return -1;
    // Parse Value
    while (gem_expect_true(IS_HEX_DIGIT(**text_line))) {
      number = (number<<4) + GET_HEX_DIGIT(**text_line);
      PARSER_NEXT_CHAR(text_line);
    }
  } else {
    // Parse sign (if any)
    bool positive = true;
    if (gem_expect_false(**text_line==PLUS)) {
      PARSER_NEXT_CHAR(text_line);
    } else if (gem_expect_false(**text_line==MINUS)) {
      positive = false;
      PARSER_NEXT_CHAR(text_line);
    }
    // Parse number
    if (gem_expect_false(!IS_DIGIT(**text_line))) return -1;
    while (gem_expect_true(IS_DIGIT(**text_line))) {
      number = (number*10) + GET_DIGIT(**text_line);
      PARSER_NEXT_CHAR(text_line);
    }
    // Add sign
    if (gem_expect_false(!positive)) number = -number;
  }
  *value = number;
  return 0;
}
GEM_INLINE error_code_t input_text_parse_double(const char** const text_line,double* const value) {
  GEM_CHECK_NULL(text_line);
  GEM_CHECK_NULL(value);
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
  if (gem_expect_false(**text_line==PLUS)) {
    PARSER_NEXT_CHAR(text_line);
  } else if (gem_expect_false(**text_line==MINUS)) {
    positive = false;
    PARSER_NEXT_CHAR(text_line);
  }
  // Parse Integer
  double integer = 0.0;
  if (gem_expect_true(IS_DIGIT(**text_line))) { // 45
    while (gem_expect_true(IS_DIGIT(**text_line))) {
      integer = (integer*10) + GET_DIGIT(**text_line);
      PARSER_NEXT_CHAR(text_line);
    }
  }
  // Dot
  if (gem_expect_true(**text_line==DOT)) { // .045
    PARSER_NEXT_CHAR(text_line);
    // Dot forces to have a decimal part
    if (gem_expect_false(!IS_DIGIT(**text_line))) return -1;
  }
  // Parse decimal
  double decimal = 0.0;
  if (gem_expect_true(IS_DIGIT(**text_line))) { // 45
    while (gem_expect_true(IS_DIGIT(**text_line))) {
      decimal = (decimal/10) + GET_DIGIT(**text_line);
      PARSER_NEXT_CHAR(text_line);
    }
  }
  // Parse exponent
  double exponent = 0.0;
  if (gem_expect_false(**text_line=='e' || **text_line=='E')) {
    PARSER_NEXT_CHAR(text_line);
    // Parse sign (if any)
    bool exponent_positive = true;
    if (gem_expect_false(**text_line==PLUS)) {
      PARSER_NEXT_CHAR(text_line);
    } else if (gem_expect_false(**text_line==MINUS)) {
      exponent_positive = false;
      PARSER_NEXT_CHAR(text_line);
    }
    // Parse number
    if (gem_expect_false(!IS_DIGIT(**text_line))) return -1;
    while (gem_expect_true(IS_DIGIT(**text_line))) {
      exponent = (exponent*10) + GET_DIGIT(**text_line);
      PARSER_NEXT_CHAR(text_line);
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
// Parsing CMD-line options
GEM_INLINE error_code_t input_text_parse_size(char* const size_text,uint64_t* const size) {
  char* text_centinel = size_text;
  double parsed_size;
  // Parse number (double/interger)
  if (input_text_parse_double((const char** const)&text_centinel,&parsed_size)) return 1;
  // Parse units
  switch (*text_centinel) {
    case 'G': /* Giga */
    case 'g':
      ++text_centinel;
      if (*text_centinel!='B' && *text_centinel!='b' && *text_centinel!='\0') return 1;
      *size = ((uint64_t)parsed_size)*1024*1024*1024;
      break;
    case 'M': /* Mega */
    case 'm':
      ++text_centinel;
      if (*text_centinel!='B' && *text_centinel!='b' && *text_centinel!='\0') return 1;
      *size = ((uint64_t)parsed_size)*1024*1024;
      break;
    case 'K': /* Kilo */
    case 'k':
      ++text_centinel;
      if (*text_centinel!='B' && *text_centinel!='b' && *text_centinel!='\0') return 1;
      *size = ((uint64_t)parsed_size)*1024;
      break;
    case 'B': /* Bytes */
      ++text_centinel;
      if (*text_centinel!='\0') return 1;
      *size = ((uint64_t)parsed_size);
      break;
    case '\0': /* Bytes */
      *size = ((uint64_t)parsed_size);
      break;
    default:
      *size = 0;
      return 1;
      break;
  }
  return 0;
}
/*
 * Tag parsing
 */
GEM_INLINE uint64_t input_text_parse_count_colons_in_field(const char* const text_line) {
  // Count number of colons in a field (delimited by TAB or SPACE)
  uint64_t i = 0, count = 0;
  while (text_line[i]!=TAB && text_line[i]!=SPACE && text_line[i]!=EOL) {
    if (text_line[i]==':') ++count;
    ++i;
  }
  return count;
}
GEM_INLINE error_code_t input_text_parse_tag(
    char** const text_line,string_t* const tag,sequence_attributes_t* const attributes) {
  GEM_CHECK_NULL(text_line);
  STRING_CHECK(tag);
  GEM_CHECK_NULL(attributes);
  // Delimit the tag
  register uint64_t i = 0;
  char* const tag_begin = *text_line;
  // Parse Tag
  PARSER_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE); // Read until first SPACE or TAB
  string_set_buffer(tag,tag_begin,(*text_line-tag_begin));
  // Add pair info and chomp /1/2/3 info (if any)
  attributes->end_info = input_text_parse_tag_chomp_pairend_info(tag);
  string_append_eos(tag);
  /*
   * Parse all extra TAG-info
   */
  while (!PARSER_IS_EOL(text_line) && **text_line!=TAB) {
    PARSER_NEXT_CHAR(text_line); // Skip space
    /*
     * CASAVA 1.8 Attributes. @EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACG
     */
    if (input_text_parse_count_colons_in_field(*text_line)==3) {
      char* const casava_tag_begin = *text_line;
      switch (casava_tag_begin[0]) {
        case '1':
          attributes->end_info = PAIRED_END1;
          break;
        case '2':
        case '3':
          attributes->end_info = PAIRED_END2;
          break;
        default:
          attributes->end_info = SINGLE_END;
          break;
      }
      PARSER_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE);
      string_set_buffer(&attributes->casava_tag,casava_tag_begin,(*text_line-casava_tag_begin));
      continue; // Next!
    }
    /*
     * Extra Tag (Consider everything else as extra information appended to the tag)
     */
    char* const extra_tag_begin = *text_line;
    PARSER_READ_UNTIL(text_line,**text_line==TAB);
    string_set_buffer(&attributes->extra_tag,extra_tag_begin,*text_line-extra_tag_begin);
    /*
     * Additional check to see if any extra attributes end in /1 /2 /3
     *     If no pair has been found yet, /1/2/3 is cut away
     *     and the tag is reset to the original tag but spaces are replaced with _
     *     E.g
     *       @SRR384920.1 HWI-ST382_0049:1:1:1217:1879/1
     *       @SRR384920.1_HWI-ST382_0049:1:1:1217:1879
     */
    const sequence_end_t end_info = input_text_parse_tag_chomp_pairend_info(&attributes->extra_tag);
    if (end_info==PAIRED_END1 || end_info==PAIRED_END2) {
      // Append to the tag an '_' plus the extra information
      string_append_char(tag,UNDERSCORE);
      string_append_string(tag,&attributes->extra_tag);
      // Replace all spaces
      const uint64_t tag_length = string_get_length(tag);
      for (i=0;i<tag_length;i++) {
        if (tag->buffer[i]==SPACE) tag->buffer[i] = UNDERSCORE;
      }
      string_append_eos(tag);
      string_clear(&attributes->extra_tag);
      // Set pair info
      attributes->end_info = end_info;
    }
  } /* while (not end of tags) */
  PARSER_NEXT_CHAR(text_line);
  return GEM_STATUS_OK;
}
/*
 * Parse the end information {/1,/2}
 */
GEM_INLINE uint64_t input_text_parse_tag_chomp_pairend_info(string_t* const tag) {
  STRING_CHECK(tag);
  const uint64_t tag_length = string_get_length(tag);
  if (tag_length>2 && *string_char_at(tag,tag_length-2)==SLASH) {
    switch (*string_char_at(tag,tag_length-1)) {
      case '1':
        string_set_length(tag,tag_length-2);
        return PAIRED_END1;
        break;
      case '2':
      case '3':
        string_set_length(tag,tag_length-2);
        return PAIRED_END2;
        break;
      default:
        break;
    }
  }
  return SINGLE_END;
}
