/*
 * PROJECT: GEMMapper
 * FILE: input_parser.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "io/input_parser.h"

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


//
///*
// * Basic InputFile-Parsing Functions
// */
//bool input_file_parse_next_char(input_file_t* const input_file) {
//  return input_file_next_char(input_file);
//}
//error_code_t input_file_parse_skip_separators(input_file_t* const input_file) {
//  if (gem_expect_false(input_file_eof(input_file))) return -1;
//  const char current_char = input_file_get_current_char(input_file);
//  if (gem_expect_false(current_char!=TAB && current_char!=SPACE)) return -1;
//  // Skip until different from TAB or SPACE
//  while (input_file_next_char(input_file)) {
//    // Read character
//    const char current_char = input_file_get_current_char(input_file);
//    if (IS_ANY_EOL(current_char) || (current_char!=TAB && current_char!=SPACE)) break;
//  }
//  return 0;
//}
//void input_file_parse_skip_chars(input_file_t* const input_file,uint64_t num_chars) {
//  while (!input_file_eof(input_file) && num_chars>0) {
//    // Read character
//    const char current_char = input_file_get_current_char(input_file);
//    if (IS_ANY_EOL(current_char)) break;
//    // Next
//    input_file_next_char(input_file);
//    --num_chars;
//  }
//}
//void input_file_parse_skip_line(input_file_t* const input_file) {
//  while (!input_file_eof(input_file)) {
//    // Read character
//    const char current_char = input_file_get_current_char(input_file);
//    if (IS_ANY_EOL(current_char)) break;
//    // Next
//    input_file_next_char(input_file);
//  }
//}
//bool input_file_parse_is_eol(input_file_t* const input_file) {
//  return input_file_eof(input_file) || IS_ANY_EOL(input_file_get_current_char(input_file));
//}
//void input_file_parse_field(input_file_t* const input_file,const char delimiter,string_t* const string) {
//  while (!input_file_eof(input_file)) {
//    // Read character
//    const char current_char = input_file_get_current_char(input_file);
//    if (IS_ANY_EOL(current_char) || current_char==delimiter) break;
//    // Append character
//    string_append_char(string,current_char);
//    // Next
//    input_file_next_char(input_file);
//  }
//  string_append_eos(string);
//}
//error_code_t input_file_parse_integer(input_file_t* const input_file,int64_t* const value) {
//  int64_t number = 0;
//  if (input_file_eof(input_file)) return -1;
//  char current_char = input_file_get_current_char(input_file);
//  // Parse sign (if any)
//  bool positive = true;
//  if (gem_expect_false(current_char==PLUS)) {
//    if (!input_file_next_char(input_file)) return -1;
//    current_char = input_file_get_current_char(input_file);
//  } else if (gem_expect_false(current_char==MINUS)) {
//    positive = false;
//    if (!input_file_next_char(input_file)) return -1;
//    current_char = input_file_get_current_char(input_file);
//  }
//  // Parse number
//  if (gem_expect_false(!IS_DIGIT(current_char))) return -1;
//  while (gem_expect_true(IS_DIGIT(current_char))) {
//    number = (number*10) + GET_DIGIT(current_char);
//    if (!input_file_next_char(input_file)) break;
//    current_char = input_file_get_current_char(input_file);
//  }
//  // Add sign
//  if (gem_expect_false(!positive)) number = -number;
//  *value = number;
//  return 0;
//}
//
///*
// * Basic Text-Parsing Functions
// */
//void input_text_parse_next_char(const char** const text_line) {
//  PARSER_NEXT_CHAR(text_line);
//}
//void input_text_parse_skip_chars(const char** const text_line,uint64_t num_chars) {
//  while ((num_chars--) > 0 && !PARSER_IS_EOL(text_line)) PARSER_NEXT_CHAR(text_line);
//}
//void input_text_parse_skip_line(const char** const text_line) {
//  PARSER_SKIP_LINE(text_line);
//}
//bool input_text_parse_is_eol(const char** const text_line) {
//  return PARSER_IS_EOL(text_line);
//}
//void input_text_parse_field(const char** const text_line,const char delimiter,string_t* const string) {
//  // Read field
//  const char* const string_begin = *text_line;
//  while (gem_expect_true(**text_line!=delimiter && !PARSER_IS_EOL(text_line))) PARSER_NEXT_CHAR(text_line);
//  // Copy string
//  if (string) string_set_buffer_const(string,string_begin,(*text_line-string_begin));
//  // Skip delimiter
//  if (**text_line==delimiter) PARSER_NEXT_CHAR(text_line);
//}
int input_text_parse_integer(const char** const text_line,int64_t* const value) {
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
int input_text_parse_double(const char** const text_line,double* const value) {
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
  double decimal = 0.0;
  if (gem_expect_true(**text_line==DOT)) { // .045
    PARSER_NEXT_CHAR(text_line);
    // Dot forces to have a decimal part
    if (gem_expect_false(!IS_DIGIT(**text_line))) return -1;
    // Parse decimal
    double decimal_div = 10.0;
    if (gem_expect_true(IS_DIGIT(**text_line))) { // 45
      while (gem_expect_true(IS_DIGIT(**text_line))) {
        decimal += (double)GET_DIGIT(**text_line)/decimal_div;
        decimal_div *= 10.0;
        PARSER_NEXT_CHAR(text_line);
      }
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
  if (exponent!=0.0) *value = *value * pow(10,exponent);
  if (!positive) *value = -(*value);
  return 0;
}
// Parsing CMD-line options
int input_text_parse_size(char* const size_text,uint64_t* const size) {
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
int input_text_parse_csv_arguments(char* const arguments,const uint64_t num_arguments,...) {
  uint64_t num_arguments_parsed = 0;
  // Start va_args
  va_list v_args;
  va_start(v_args,num_arguments);
  // Start parsing
  char *opt = strtok(arguments,",");
  while (opt!=NULL && num_arguments_parsed<num_arguments) {
    char** const arg = va_arg(v_args,char**);
    *arg = opt;
    opt = strtok(NULL,",");
    ++num_arguments_parsed;
  }
  // End va_args
  va_end(v_args);
  return num_arguments_parsed;
}
int input_text_parse_extended_uint64(char* const argument,uint64_t* const value) {
  // Textual
  if (gem_strcaseeq(argument,"all")) { *value = UINT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"inf")) { *value = UINT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"infinite"))  { *value = UINT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"unlimited")) { *value = UINT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"none")) { *value = 0; return 0; }
  if (gem_strcaseeq(argument,"zero")) { *value = 0; return 0; }
  if (gem_strcaseeq(argument,"null")) { *value = 0; return 0; }
  // Number
  return input_text_parse_integer((const char** const)&argument,(int64_t*)value);
}
int input_text_parse_extended_int64(char* const argument,int64_t* const value) {
  // Textual
  if (gem_strcaseeq(argument,"all")) { *value = INT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"inf")) { *value = INT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"infinite")) { *value = INT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"-inf")) { *value = INT64_MIN; return 0; }
  if (gem_strcaseeq(argument,"-infinite")) { *value = INT64_MIN; return 0; }
  if (gem_strcaseeq(argument,"none")) { *value = 0; return 0; }
  if (gem_strcaseeq(argument,"zero")) { *value = 0; return 0; }
  // Number
  return input_text_parse_integer((const char** const)&argument,value);
}
int input_text_parse_extended_double(char* const argument,double* const value) {
  // Textual (use int64_t limits)
  if (gem_strcaseeq(argument,"all")) { *value = INT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"inf")) { *value = INT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"infinite")) { *value = INT64_MAX; return 0; }
  if (gem_strcaseeq(argument,"-inf")) { *value = INT64_MIN; return 0; }
  if (gem_strcaseeq(argument,"-infinite")) { *value = INT64_MIN; return 0; }
  if (gem_strcaseeq(argument,"none")) { *value = 0; return 0; }
  if (gem_strcaseeq(argument,"zero")) { *value = 0; return 0; }
  // Number
  return input_text_parse_double((const char** const)&argument,value);
}
bool input_text_parse_extended_bool(char* const argument) {
  if (argument==NULL) {
    return true;
  } else {
    if (gem_streq(argument,"true") ||
        gem_streq(argument,"yes")  ||
        gem_streq(argument,"allow") ||
        gem_streq(argument,"enabled")) {
      return true;
    } else {
      return false;
    }
  }
}
uint64_t input_text_parse_count_colons_in_field(const char* const text_line) {
  // Count number of colons in a field (delimited by TAB or SPACE)
  uint64_t i = 0, count = 0;
  while (text_line[i]!=TAB && text_line[i]!=SPACE && text_line[i]!=EOL) {
    if (text_line[i]==':') ++count;
    ++i;
  }
  return count;
}

