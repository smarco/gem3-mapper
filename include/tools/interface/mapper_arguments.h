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

#ifndef MAPPER_ARGUMENTS_H_
#define MAPPER_ARGUMENTS_H_

#include "utils/essentials.h"
#include "utils/options_menu.h"
#include "mapper/mapper_parameters.h"

/*
 * Mapper options Menu
 */
extern option_t gem_mapper_options[];
extern char* gem_mapper_groups[];

/*
 * Mapper Error
 */
#define mapper_error_msg(error_msg,args...) \
  fprintf(stderr,"GEM-Mapper error:\n> "error_msg"\n",##args); \
  exit(1)
#define mapper_cond_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      mapper_error_msg(error_msg,##args); \
      exit(1); \
    } \
  } while (0)

/*
 * Mapper Usage
 */
void gem_mapper_print_usage(const option_visibility_t visibility_level);

/*
 * Mapper Arguments Parsing
 */
void gem_mapper_parse_arguments(
    int argc,
    char** argv,
    mapper_parameters_t* const parameters,
    char* const gem_version);

#endif /* MAPPER_ARGUMENTS_H_ */
