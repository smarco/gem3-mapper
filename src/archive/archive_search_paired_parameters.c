/*
 * PROJECT: GEMMapper
 * FILE: archive_search_paired_parameters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "archive/archive_search_paired_parameters.h"

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
  search_paired_parameters->min_initial_template_estimation = UINT64_MAX;
  search_paired_parameters->max_initial_template_estimation = UINT64_MAX;
}
