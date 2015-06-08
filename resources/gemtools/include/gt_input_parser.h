/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_parser.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Thasso Griebel <thasso.griebel@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_PARSER_H_
#define GT_INPUT_PARSER_H_

#include "gt_essentials.h"
#include "gt_attributes.h"

/*
 * Parsing error/state codes
 */
#define GT_PE_BAD_CHARACTER 10

#define GT_PE_BAD_TRIM_READ_STRING_LENGTH 30
#define GT_PE_BAD_TRIM_QUAL_STRING_LENGTH 31

#define GT_INPUT_PARSER_IS_SAM_ATTRIBUTE(text_line) \
  (!gt_is_end_of_field(**text_line) && \
   !gt_is_end_of_field(*(*text_line+1)) && \
   (*(*text_line+2))==COLON && \
   !gt_is_end_of_field(*(*text_line+3)) && \
   (*(*text_line+4))==COLON)

/*
 * Basic Parsing Functions
 */
GT_INLINE void gt_input_parse_next_char(const char** const text_line);
GT_INLINE void gt_input_parse_skip_chars(const char** const text_line,uint64_t num_chars);
GT_INLINE gt_status gt_input_parse_eol(const char** const text_line);
GT_INLINE void gt_input_parse_field(const char** const text_line,const char delimiter,gt_string* const string);
GT_INLINE gt_status gt_input_parse_integer(const char** const text_line,int64_t* const value);
GT_INLINE gt_status gt_input_parse_double(const char** const text_line,double* const value);

/*
 * Parse any SAM-like attribute
 */
GT_INLINE gt_status gt_input_parse_sam_optional_field(const char** const text_line,gt_attributes* const attributes);

/*
 * Generic Tag Parser & Utils
 */
GT_INLINE gt_status gt_input_parse_tag(const char** const text_line,gt_string* const tag,gt_attributes* const attributes);
GT_INLINE uint64_t gt_input_parse_tag_chomp_pairend_info(gt_string* const tag);

/*
 * Internal Building Blocks for parsing
 */
#define GT_NEXT_CHAR(text_line) ++(*text_line)
#define GT_IS_EOL(text_line) gt_expect_false((**text_line)==EOL || (**text_line)==EOS)
#define GT_READ_UNTIL(text_line,test) \
  while (gt_expect_true(!(test) && !GT_IS_EOL(text_line))) { \
    GT_NEXT_CHAR(text_line); \
  }
#define GT_PARSE_HEX_OR_DEC(text_line,number) \
  number=0; \
  if (**text_line=='0' && (*(*text_line+1)=='x' || *(*text_line+1)=='X')) { \
    *text_line+=2; \
    while (gt_expect_true(gt_is_hex_digit(**text_line))) { \
      number = (number<<4) + gt_get_hex_cipher(**text_line); \
      GT_NEXT_CHAR(text_line); \
    } \
  } else { \
    while (gt_expect_true(gt_is_number(**text_line))) { \
      number = (number*10) + gt_get_cipher(**text_line); \
      GT_NEXT_CHAR(text_line); \
    } \
  }
#define GT_PARSE_NUMBER(text_line,number) \
  number = 0; \
  while (gt_expect_true(gt_is_number(**text_line))) { \
    number = (number*10) + gt_get_cipher(**text_line); \
    GT_NEXT_CHAR(text_line); \
  }
#define GT_PARSE_SIGNED_NUMBER_BLOCK(text_line,number) { \
  bool is_negative; \
  switch ((**text_line)) { \
    case PLUS: is_negative = false; break; \
    case MINUS: is_negative = true; break; \
    default: is_negative = false; break; \
  }
#define GT_PARSE_SIGNED_NUMBER_END_BLOCK(number) \
  if (is_negative) number = -number; \
}
#define GT_SKIP_LINE(text_line) \
  while (!GT_IS_EOL(text_line)) { \
    ++(*text_line); \
  }

#endif /* GT_INPUT_PARSER_H_ */
