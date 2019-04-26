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
 */

#ifndef OPTIONS_MENU_H_
#define OPTIONS_MENU_H_

#include "utils/essentials.h"
#include <getopt.h>

/*
 * Options (Tools Menu)
 */
typedef enum { NO_ARGUMENT=no_argument, REQUIRED=required_argument, OPTIONAL=optional_argument } option_type;
typedef enum { TYPE_NONE, TYPE_INT, TYPE_FLOAT, TYPE_CHAR, TYPE_STRING, TYPE_BOOL } option_data_t;
typedef enum { VISIBILITY_USER=0, VISIBILITY_ADVANCED=1, VISIBILITY_DEVELOPER=2 } option_visibility_t;
typedef struct {
  int option_id;                         // Integer ID or short character option
  char* long_option;                     // Long string option
  option_type option_type;               // Argument type
  option_data_t argument_type;           // Argument Data type
  uint64_t group_id;                     // Label of the group it belongs to (zero if none)
  option_visibility_t option_visibility; // Enable/Disable option
  char* command_info;                    // Extra command line syntax info
  char* description;                     // Brief description
} option_t;

uint64_t options_get_num_options(const option_t* const options_menu);
struct option* options_adaptor_getopt(const option_t* const options_menu);
string_t* options_adaptor_getopt_short(const option_t* const options_menu);
void options_fprint_menu(
    FILE* const stream,
    const option_t* const options_menu,
    char* groups_menu[],
    const bool print_description,
    const option_visibility_t visibility_level);

#define GEM_OPTIONS_ITERATE_BEGIN(tool_options,option_id) \
  struct option* __getopt = gt_options_adaptor_getopt(tool_options); \
  gt_string* const __short_getopt = gt_options_adaptor_getopt_short(tool_options); \
  int option_id, option_index; \
  while (true) { \
    /* Get option &  Select case */ \
    if ((option_id=getopt_long(argc,argv,gt_string_get_string(__short_getopt),__getopt,&option_index))==-1) break; \
    switch (option_id)

#define GEM_OPTIONS_ITERATE_END \
    string_delete(__short_getopt)

/*
 * Useful option parsing tools
 */
#define options_parse_bool(str_argument) (gem_streq(str_argument,"true") || gem_streq(str_argument,"yes"))

#endif /* OPTIONS_MENU_H_ */

