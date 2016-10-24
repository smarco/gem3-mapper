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

#include "archive/search/archive_search_pe_parameters.h"

/*
 * Paired-Search Parameters
 */
void search_paired_parameters_init(search_paired_parameters_t* const search_paired_parameters) {
  search_paired_parameters->paired_end_search = false;
  search_paired_parameters->paired_end_extension_shortcut = false;
  search_paired_parameters->paired_end_extension_recovery = true;
  /* Concordant Orientation */
  search_paired_parameters->pair_discordant_search = pair_discordant_search_only_if_no_concordant;
  search_paired_parameters->pair_orientation[pair_orientation_FR] = pair_relation_concordant;
  search_paired_parameters->pair_orientation[pair_orientation_RF] = pair_relation_discordant;
  search_paired_parameters->pair_orientation[pair_orientation_FF] = pair_relation_invalid;
  search_paired_parameters->pair_orientation[pair_orientation_RR] = pair_relation_invalid;
  /* Pair allowed lay-outs */
  search_paired_parameters->pair_layout[pair_layout_separate] = pair_relation_concordant;
  search_paired_parameters->pair_layout[pair_layout_overlap] = pair_relation_concordant;
  search_paired_parameters->pair_layout[pair_layout_contain] = pair_relation_discordant;
  /* Template length constraints */
  search_paired_parameters->min_template_length = 0;
  search_paired_parameters->max_template_length = 10000;
  search_paired_parameters->template_length_estimation_min = 0;
  search_paired_parameters->template_length_estimation_max = 1000;
  search_paired_parameters->template_length_estimation_samples = 100;
}
