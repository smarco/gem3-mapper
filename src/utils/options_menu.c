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

#include "utils/options_menu.h"

uint64_t options_get_num_options(const option_t* const options_menu) {
  uint64_t num_options = 0, i = 0;
  while (options_menu[i++].option_id != 0) ++num_options;
  return num_options;
}
struct option* options_adaptor_getopt(const option_t* const options_menu) {
  const uint64_t num_options = options_get_num_options(options_menu);
  struct option* getopt_menu = mm_malloc(sizeof(struct option)*(num_options+1));
  // Adapt all the records
  uint64_t i = 0;
  for (i=0;i<num_options;++i) {
    getopt_menu[i].name = options_menu[i].long_option;
    getopt_menu[i].has_arg = options_menu[i].option_type;
    getopt_menu[i].flag = 0;
    getopt_menu[i].val = options_menu[i].option_id;
  }
  getopt_menu[num_options].name = 0;
  getopt_menu[num_options].has_arg = 0;
  getopt_menu[num_options].flag = 0;
  getopt_menu[num_options].val = 0;
  return getopt_menu;
}
string_t* options_adaptor_getopt_short(const option_t* const options_menu) {
  const uint64_t num_options = options_get_num_options(options_menu);
  string_t* const getopt_menu_short = mm_alloc(string_t);
  string_init(getopt_menu_short,2*num_options,NULL);
  // Adapt all the short options
  uint64_t i = 0;
  for (i=0;i<num_options;++i) {
    const char short_option = options_menu[i].option_id;
    if (options_menu[i].option_id<128 && IS_ALPHANUMERIC(short_option)) {
      string_append_char(getopt_menu_short,short_option);
      if (options_menu[i].option_type==REQUIRED) {
        string_append_char(getopt_menu_short,COLON);
      } else if (options_menu[i].option_type==OPTIONAL) {
        string_append_char(getopt_menu_short,COLON);
        string_append_char(getopt_menu_short,COLON);
      }
    }
  }
  string_append_eos(getopt_menu_short);
  return getopt_menu_short;
}
void options_fprint_menu(
    FILE* const stream,
    const option_t* const options_menu,
    char* groups_menu[],
    const bool print_description,
    const option_visibility_t visibility_level) {
  const uint64_t num_options = options_get_num_options(options_menu);
  int64_t i, last_group = -1;
  for (i=0;i<num_options;++i) {
    if (options_menu[i].option_visibility > visibility_level) continue;
    // Print group (if not printed yet)
    if (last_group!=options_menu[i].group_id) {
      fprintf(stream,"    [%s]\n",groups_menu[options_menu[i].group_id]);
      last_group=options_menu[i].group_id;
    }
    // Print Long Option
    fprintf(stream,"      --%s",options_menu[i].long_option);
    // Print Short Option (if it has)
    const char short_option = options_menu[i].option_id;
    if (options_menu[i].option_id<128 && IS_ALPHANUMERIC(short_option)) {
      fprintf(stream,"|-%c",short_option);
    }
    // Print extra command line syntax info
    fprintf(stream," %s",options_menu[i].command_info);
    // Print description (@print_description)
    if (print_description && !gem_streq(options_menu[i].description,"")) {
      fprintf(stream," %s\n",options_menu[i].description);
    } else {
      fprintf(stream,"\n");
    }
  }
}
