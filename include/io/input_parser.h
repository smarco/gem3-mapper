/*
 * PROJECT: GEMMapper
 * FILE: input_parser.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef INPUT_PARSER_H_
#define INPUT_PARSER_H_

#include "utils/essentials.h"
#include "io/input_file.h"
#include "data_structures/sequence.h"

/*
 * Errors
 */
#define GEM_ERROR_PARSING_SIZE "Error parsing %s. '%s' not a valid size (Eg. 2GB)"

///*
// * Basic Input File Parsing Functions
// */
//bool input_file_parse_next_char(input_file_t* const restrict input_file);
//error_code_t input_file_parse_skip_separators(input_file_t* const restrict input_file);
//void input_file_parse_skip_chars(input_file_t* const restrict input_file,uint64_t num_chars);
//void input_file_parse_skip_line(input_file_t* const restrict input_file);
//bool input_file_parse_is_eol(input_file_t* const restrict input_file);
//void input_file_parse_field(input_file_t* const restrict input_file,const char delimiter,string_t* const restrict string);
//error_code_t input_file_parse_integer(input_file_t* const restrict input_file,int64_t* const restrict value);
//error_code_t input_file_parse_double(input_file_t* const restrict input_file,double* const restrict value);
//
/*
 * Basic Text Parsing Functions
 */
//void input_text_parse_next_char(const char** const restrict text_line);
//void input_text_parse_skip_chars(const char** const restrict text_line,uint64_t num_chars);
//void input_text_parse_skip_line(const char** const restrict text_line);
//bool input_text_parse_is_eol(const char** const restrict text_line);
//void input_text_parse_field(const char** const restrict text_line,const char delimiter,string_t* const restrict string);
int input_text_parse_integer(const char** const restrict text_line,int64_t* const restrict value);
int input_text_parse_double(const char** const restrict text_line,double* const restrict value);
int input_text_parse_size(char* const restrict size_text,uint64_t* const restrict size);
int input_text_parse_csv_arguments(char* const restrict arguments,const uint64_t num_arguments,...);
int input_text_parse_extended_uint64(char* const restrict argument,uint64_t* const restrict value);
int input_text_parse_extended_int64(char* const restrict argument,int64_t* const restrict value);
int input_text_parse_extended_double(char* const restrict argument,double* const restrict value);
bool input_text_parse_extended_bool(char* const restrict argument);
uint64_t input_text_parse_count_colons_in_field(const char* const restrict text_line);

#endif /* INPUT_PARSER_H_ */
