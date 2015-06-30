/*
 * PROJECT: GEMMapper
 * FILE: paired_search_parameters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PAIRED_SEARCH_PARAMETERS_H_
#define PAIRED_SEARCH_PARAMETERS_H_

#include "essentials.h"

typedef enum {
  pair_discordant_search_always,
  pair_discordant_search_only_if_no_concordant,
  pair_discordant_search_never
} pair_discordant_search_t;
typedef enum {
  pair_layout_separate,
  pair_layout_overlap,
  pair_layout_contain,
  pair_layout_dovetail,
} pair_layout_t;
typedef enum {
  pair_orientation_concordant,
  pair_orientation_discordant,
  pair_orientation_invalid
} pair_orientation_t;
typedef struct {
  pair_discordant_search_t pair_discordant_search;
  uint64_t max_paired_matches;
  /* Pair Orientation */
  pair_orientation_t pair_orientation_FR;
  pair_orientation_t pair_orientation_RF;
  pair_orientation_t pair_orientation_FF;
  pair_orientation_t pair_orientation_RR;
  /* Pair allowed lay-outs */
  bool pair_layout_separate;
  bool pair_layout_overlap;
  bool pair_layout_contain;
  bool pair_layout_dovetail;
  /* Template length constraints */
  uint64_t min_template_length;
  uint64_t max_template_length;
} paired_search_parameters_t;

/*
 * Paired-Search Parameters
 */
GEM_INLINE void paired_search_parameters_init(paired_search_parameters_t* const paired_search_parameters);

#endif /* PAIRED_SEARCH_PARAMETERS_H_ */
