/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Input module provides simple parsing functions (i.e. fields, numbers, ...)
 */

#ifndef INPUT_PARSER_H_
#define INPUT_PARSER_H_

#include "utils/essentials.h"
#include "io/input_file.h"
#include "text/sequence.h"

/*
 * Errors
 */
#define GEM_ERROR_PARSING_SIZE "Error parsing %s. '%s' not a valid size (Eg. 2GB)"

///*
// * Basic Input File Parsing Functions
// */
//bool input_file_parse_next_char(input_file_t* const input_file);
//error_code_t input_file_parse_skip_separators(input_file_t* const input_file);
//void input_file_parse_skip_chars(input_file_t* const input_file,uint64_t num_chars);
//void input_file_parse_skip_line(input_file_t* const input_file);
//bool input_file_parse_is_eol(input_file_t* const input_file);
//void input_file_parse_field(input_file_t* const input_file,const char delimiter,string_t* const string);
//error_code_t input_file_parse_integer(input_file_t* const input_file,int64_t* const value);
//error_code_t input_file_parse_double(input_file_t* const input_file,double* const value);
//
/*
 * Basic Text Parsing Functions
 */
//void input_text_parse_next_char(const char** const text_line);
//void input_text_parse_skip_chars(const char** const text_line,uint64_t num_chars);
//void input_text_parse_skip_line(const char** const text_line);
//bool input_text_parse_is_eol(const char** const text_line);
//void input_text_parse_field(const char** const text_line,const char delimiter,string_t* const string);
int input_text_parse_integer(const char** const text_line,int64_t* const value);
int input_text_parse_double(const char** const text_line,double* const value);
int input_text_parse_size(char* const size_text,uint64_t* const size);
int input_text_parse_csv_arguments(char* const arguments,const uint64_t num_arguments,...);
int input_text_parse_extended_uint64(char* const argument,uint64_t* const value);
int input_text_parse_extended_int64(char* const argument,int64_t* const value);
int input_text_parse_extended_double(char* const argument,double* const value);
bool input_text_parse_extended_bool(char* const argument);
uint64_t input_text_parse_count_colons_in_field(const char* const text_line);

#endif /* INPUT_PARSER_H_ */
