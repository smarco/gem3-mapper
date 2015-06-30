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
GEM_INLINE void paired_search_parameters_init(paired_search_parameters_t* const paired_search_parameters) {
  paired_search_parameters->max_paired_matches = 1000;
  /* Concordant Orientation */
  paired_search_parameters->pair_discordant_search = pair_discordant_search_only_if_no_concordant;
  paired_search_parameters->pair_orientation_FR = pair_orientation_concordant;
  paired_search_parameters->pair_orientation_RF = pair_orientation_invalid;
  paired_search_parameters->pair_orientation_FF = pair_orientation_invalid;
  paired_search_parameters->pair_orientation_RR = pair_orientation_invalid;
  /* Pair allowed lay-outs */
  paired_search_parameters->pair_layout_separate = true;
  paired_search_parameters->pair_layout_overlap = true;
  paired_search_parameters->pair_layout_contain = true;
  paired_search_parameters->pair_layout_dovetail = true;
  /* Template length constraints */
  paired_search_parameters->min_template_length = 0;
  paired_search_parameters->max_template_length = 10000;
}
