/*
 * PROJECT: GEMMapper
 * FILE: matches_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches/matches_metrics.h"

/*
 * Setup
 */
void matches_metrics_init(matches_metrics_t* const metrics) {
  // Aggregated
  metrics->total_matches_sampled = 0;
  metrics->accepted_candidates = 0;
  // Minimums
  metrics->min1_counter_value = UINT32_MAX;
  metrics->min2_counter_value = UINT32_MAX;
  metrics->min1_edit_distance = UINT32_MAX;
  metrics->min2_edit_distance = UINT32_MAX;
  metrics->max1_swg_score = INT32_MIN;
  metrics->max2_swg_score = INT32_MIN;
  // Template length
  metrics->min1_template_length_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
  metrics->min2_template_length_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
  // MAPQ Score
  metrics->mapq = 0;
}
/*
 * Accessors
 */
uint64_t matches_metrics_get_min_distance(matches_metrics_t* const metrics) {
  return metrics->min1_counter_value;
}
uint64_t matches_metrics_get_min_edit_distance(matches_metrics_t* const metrics) {
  return metrics->min1_edit_distance;
}
int32_t matches_metrics_get_max_swg_score(matches_metrics_t* const metrics) {
  return metrics->max1_swg_score;
}
void matches_metrics_add_accepted_candidates(matches_metrics_t* const metrics,const uint64_t num_candidates) {
  metrics->accepted_candidates += num_candidates;
}
void matches_metrics_set_mapq(matches_metrics_t* const metrics,const uint8_t mapq) {
  metrics->mapq = mapq;
}
void matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t distance,
    const uint64_t edit_distance,
    const int32_t swg_score) {
  // Samples
  ++(matches_metrics->total_matches_sampled);
  // Min counter
  if (distance <= matches_metrics->min1_counter_value) {
    matches_metrics->min2_counter_value = matches_metrics->min1_counter_value;
    matches_metrics->min1_counter_value = distance;
  } else if (distance < matches_metrics->min2_counter_value) {
    matches_metrics->min2_counter_value = distance;
  }
  // Min Edit-distance
  if (edit_distance <= matches_metrics->min1_edit_distance) {
    matches_metrics->min2_edit_distance = matches_metrics->min1_edit_distance;
    matches_metrics->min1_edit_distance = edit_distance;
  } else if (edit_distance < matches_metrics->min2_edit_distance) {
    matches_metrics->min2_edit_distance = edit_distance;
  }
  // Max SGW-score
  if (swg_score >= matches_metrics->max1_swg_score) {
    matches_metrics->max2_swg_score = matches_metrics->max1_swg_score;
    matches_metrics->max1_swg_score = swg_score;
  } else if (swg_score > matches_metrics->max2_swg_score) {
    matches_metrics->max2_swg_score = swg_score;
  }
}
void paired_matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t distance,
    const uint64_t edit_distance,
    const int32_t swg_score,
    const double template_length_sigma) {
  // Update metrics
  matches_metrics_update(matches_metrics,distance,edit_distance,swg_score);
  // Min TemplateSize-Sigma
  if (template_length_sigma <= matches_metrics->min1_template_length_sigma) {
    matches_metrics->min2_template_length_sigma = matches_metrics->min1_template_length_sigma;
    matches_metrics->min1_template_length_sigma = template_length_sigma;
  } else if (template_length_sigma < matches_metrics->min2_template_length_sigma) {
    matches_metrics->min2_template_length_sigma = template_length_sigma;
  }
}
/*
 * Display
 */
void matches_metrics_print(FILE* const stream,matches_metrics_t* const matches_metrics) {
  tab_fprintf(stream,"[GEM]>Metrics\n");
  tab_fprintf(stream,"  => Aggregated \n");
  tab_fprintf(stream,"    => Total.matches.sampled %lu\n",matches_metrics->total_matches_sampled);
  tab_fprintf(stream,"    => Total.accepted.matches %lu\n",matches_metrics->accepted_candidates);
  tab_fprintf(stream,"  => Minimums\n");
  tab_fprintf(stream,"    => Min1.counter.value %lu\n",matches_metrics->min1_counter_value);
  tab_fprintf(stream,"    => Min2.counter.value %lu\n",matches_metrics->min2_counter_value);
  tab_fprintf(stream,"    => Min1.edit.distance %lu\n",matches_metrics->min1_edit_distance);
  tab_fprintf(stream,"    => Min2.edit.distance %lu\n",matches_metrics->min2_edit_distance);
  tab_fprintf(stream,"    => Max1.swg.score %ld\n",matches_metrics->max1_swg_score);
  tab_fprintf(stream,"    => Max2.swg.score %ld\n",matches_metrics->max2_swg_score);
  tab_fprintf(stream,"  => PE.Specific\n");
  tab_fprintf(stream,"    => Min1.template.sigma %2.3f\n",matches_metrics->min1_template_length_sigma);
  tab_fprintf(stream,"    => Min2.template.sigma %2.3f\n",matches_metrics->min2_template_length_sigma);
  tab_fprintf(stream,"  => MAPQ\n");
  tab_fprintf(stream,"    => MAPQ %d\n",matches_metrics->mapq);
}
