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
  pair_discordant_search_t pair_discordant_search;
  uint64_t max_paired_matches;
  /* Pair Relation Allowed {concordant,discordant,invalid} */
  pair_relation_t pair_orientation[4];
  /* Pair allowed lay-outs */
  pair_relation_t pair_layout[3];
  /* Template length constraints */
  uint64_t min_template_length;
  uint64_t max_template_length;
} paired_search_parameters_t;

/*
 * Paired-Search Parameters
 */
void paired_search_parameters_init(paired_search_parameters_t* const paired_search_parameters);

#endif /* PAIRED_SEARCH_PARAMETERS_H_ */