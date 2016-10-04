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
 *   Archive-Search Paired-End Parameters (initializers & handlers)
 */

#ifndef ARCHIVE_SEARCH_PAIRED_PARAMETERS_H_
#define ARCHIVE_SEARCH_PAIRED_PARAMETERS_H_

#include "utils/essentials.h"

typedef enum {
  pair_discordant_search_always,
  pair_discordant_search_only_if_no_concordant,
  pair_discordant_search_never
} pair_discordant_search_t;
typedef enum {
  pair_layout_separate = 0,
  pair_layout_overlap  = 1,
  pair_layout_contain  = 2,
} pair_layout_t;
typedef enum {
  pair_orientation_FR = 0,
  pair_orientation_RF = 1,
  pair_orientation_FF = 2,
  pair_orientation_RR = 3,
} pair_orientation_t;
typedef enum {
  pair_relation_concordant,
  pair_relation_discordant,
  pair_relation_invalid
} pair_relation_t;
typedef struct {
  bool paired_end_search;
  bool paired_end_extension_shortcut;
  bool paired_end_extension_recovery;
  pair_discordant_search_t pair_discordant_search;
  /* Pair Relation Allowed {concordant,discordant,invalid} */
  pair_relation_t pair_orientation[4];
  /* Pair allowed lay-outs */
  pair_relation_t pair_layout[3];
  /* Template length constraints */
  uint64_t min_template_length;
  uint64_t max_template_length;
  uint64_t template_length_estimation_min;
  uint64_t template_length_estimation_max;
  uint64_t template_length_estimation_samples;
} search_paired_parameters_t;

/*
 * Paired-Search Parameters
 */
void search_paired_parameters_init(search_paired_parameters_t* const search_paired_parameters);

#endif /* ARCHIVE_SEARCH_PAIRED_PARAMETERS_H_ */
