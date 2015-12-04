/*
 * PROJECT: GEMMapper
 * FILE: paired_search_parameters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "paired_search_parameters.h"

/*
 * Paired-Search Parameters
 */
void paired_search_parameters_init(paired_search_parameters_t* const paired_search_parameters) {
  paired_search_parameters->paired_end_search = false;
  /* Concordant Orientation */
  paired_search_parameters->pair_discordant_search = pair_discordant_search_only_if_no_concordant;
  paired_search_parameters->pair_orientation[pair_orientation_FR] = pair_relation_concordant;
  paired_search_parameters->pair_orientation[pair_orientation_RF] = pair_relation_discordant;
  paired_search_parameters->pair_orientation[pair_orientation_FF] = pair_relation_invalid;
  paired_search_parameters->pair_orientation[pair_orientation_RR] = pair_relation_invalid;
  /* Pair allowed lay-outs */
  paired_search_parameters->pair_layout[pair_layout_separate] = pair_relation_concordant;
  paired_search_parameters->pair_layout[pair_layout_overlap] = pair_relation_concordant;
  paired_search_parameters->pair_layout[pair_layout_contain] = pair_relation_discordant;
  /* Template length constraints */
  paired_search_parameters->min_template_length = 0;
  paired_search_parameters->max_template_length = 10000;
}
