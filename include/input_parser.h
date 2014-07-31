/*
 * PROJECT: GEMMapper
 * FILE: input_parser.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef INPUT_PARSER_H_
#define INPUT_PARSER_H_

#include "essentials.h"
#include "sequence.h"
#include "input_file.h"

/*
 * Checkers
 */
#define TEXT_LINE_CHECK(text_line) \
  GEM_CHECK_NULL(text_line); \
  GEM_CHECK_NULL(*text_line)

/*
 * Basic Input File Parsing Functions
 */
GEM_INLINE bool input_file_parse_next_char(input_file_t* const input_file);
GEM_INLINE error_code_t input_file_parse_skip_separators(input_file_t* const input_file);
GEM_INLINE void input_file_parse_skip_chars(input_file_t* const input_file,uint64_t num_chars);
GEM_INLINE void input_file_parse_skip_line(input_file_t* const input_file);
GEM_INLINE bool input_file_parse_is_eol(input_file_t* const input_file);
GEM_INLINE void input_file_parse_field(input_file_t* const input_file,const char delimiter,string_t* const string);
GEM_INLINE error_code_t input_file_parse_integer(input_file_t* const input_file,int64_t* const value);
GEM_INLINE error_code_t input_file_parse_double(input_file_t* const input_file,double* const value);

/*
 * Basic Text Parsing Functions
 */
GEM_INLINE void input_text_parse_next_char(const char** const text_line);
GEM_INLINE void input_text_parse_skip_chars(const char** const text_line,uint64_t num_chars);
GEM_INLINE void input_text_parse_skip_line(const char** const text_line);
GEM_INLINE bool input_text_parse_is_eol(const char** const text_line);
GEM_INLINE void input_text_parse_field(const char** const text_line,const char delimiter,string_t* const string);
GEM_INLINE error_code_t input_text_parse_integer(const char** const text_line,int64_t* const value);
GEM_INLINE error_code_t input_text_parse_double(const char** const text_line,double* const value);
GEM_INLINE error_code_t input_text_parse_size(char* const size_text,uint64_t* const size);

/*
 * Tag Parser
 */
GEM_INLINE error_code_t input_text_parse_tag(
    char** const text_line,string_t* const tag,sequence_attributes_t* const attributes);
GEM_INLINE uint64_t input_text_parse_tag_chomp_pairend_info(string_t* const tag);

/*
 * Building Blocks for parsing
 */
#define PARSER_IS_EOL(text_line) gem_expect_false((**text_line)==EOL)
#define PARSER_NEXT_CHAR(text_line) ++(*text_line)
#define PARSER_READ_UNTIL(text_line,test) \
  while (gem_expect_true(!(test) && !PARSER_IS_EOL(text_line))) { \
    PARSER_NEXT_CHAR(text_line); \
  }
#define PARSER_SKIP_LINE(text_line) \
  while (!PARSER_IS_EOL(text_line)) { \
    ++(*text_line); \
  }

/*
 * Errors
 */
#define GEM_ERROR_PARSING_SIZE "Error parsing %s. '%s' not a valid size (Eg. 2GB)"

#endif /* INPUT_PARSER_H_ */
