/*
 * PROJECT: GEMMapper
 * FILE: matches_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches_metrics.h"

/*
 * Setup
 */
void matches_metrics_init(matches_metrics_t* const metrics) {
  // Matches metrics
  metrics->total_matches_sampled = 0;
  metrics->min1_counter_value = UINT32_MAX;
  metrics->min2_counter_value = UINT32_MAX;
  metrics->max_counter_value = 0;
  metrics->min1_edit_distance = UINT32_MAX;
  metrics->min2_edit_distance = UINT32_MAX;
  metrics->max1_swg_score = INT32_MIN;
  metrics->max2_swg_score = INT32_MIN;
  metrics->subdominant_candidates = 0;
  metrics->min1_template_length_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
  metrics->min2_template_length_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
}
/*
 * Accessors
 */
uint64_t matches_metrics_get_min_distance(matches_metrics_t* const metrics) {
  return metrics->min1_counter_value;
}
uint64_t matches_metrics_get_max_distance(matches_metrics_t* const metrics) {
  return metrics->max_counter_value;
}
uint64_t matches_metrics_get_min_edit_distance(matches_metrics_t* const metrics) {
  return metrics->min1_edit_distance;
}
int32_t matches_metrics_get_max_swg_score(matches_metrics_t* const metrics) {
  return metrics->max1_swg_score;
}
/*
 * Update
 */
void matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t distance,const uint64_t edit_distance,const int32_t swg_score) {
  // Samples
  ++(matches_metrics->total_matches_sampled);
  // Max counter
  matches_metrics->max_counter_value = MAX(matches_metrics->max_counter_value,distance);
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
void matches_metrics_pe_update(
    matches_metrics_t* const matches_metrics,const uint64_t distance,
    const uint64_t edit_distance,const int32_t swg_score,
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
void matches_metrics_inc_subdominant_candidates(matches_metrics_t* const metrics) {
  ++(metrics->subdominant_candidates);
}
void matches_metrics_dec_subdominant_candidates(matches_metrics_t* const metrics) {
  --(metrics->subdominant_candidates);
}
/*
 * Display
 */
void matches_metrics_print(
    FILE* const stream,matches_metrics_t* const matches_metrics,
    const bool print_paired_metrics) {
  tab_fprintf(stream,"[GEM]>Metrics\n");
  tab_fprintf(stream,"  => Total.matches.sampled %lu\n",matches_metrics->total_matches_sampled);
  tab_fprintf(stream,"  => Total.subdominant.matches %lu\n",matches_metrics->subdominant_candidates);
  tab_fprintf(stream,"  => Min1.counter.value %lu\n",matches_metrics->min1_counter_value);
  tab_fprintf(stream,"  => Min2.counter.value %lu\n",matches_metrics->min2_counter_value);
  tab_fprintf(stream,"  => Max.counter.value %lu\n",matches_metrics->max_counter_value);
  tab_fprintf(stream,"  => Min1.edit.distance %lu\n",matches_metrics->min1_edit_distance);
  tab_fprintf(stream,"  => Min2.edit.distance %lu\n",matches_metrics->min2_edit_distance);
  tab_fprintf(stream,"  => Max1.swg.score %ld\n",matches_metrics->max1_swg_score);
  tab_fprintf(stream,"  => Max2.swg.score %ld\n",matches_metrics->max2_swg_score);
  if (print_paired_metrics) {
  tab_fprintf(stream,"  => Min1.template.sigma %f\n",matches_metrics->min1_template_length_sigma);
  tab_fprintf(stream,"  => Min2.template.sigma %f\n",matches_metrics->min2_template_length_sigma);
  }
}
void matches_predictors_print_basic_fields(
    FILE* const stream,matches_predictors_t* const predictors,
    const char* const tag,const uint8_t mapq_score) {
  // tag
  fprintf(stream,"%s\t",tag);
  // mapq_score
  fprintf(stream,"%d\t",mapq_score);
  // edit
  fprintf(stream,"%f\t",predictors->first_map_edit_distance_norm);
  // sub_edit
  fprintf(stream,"%f\t",predictors->subdominant_edit_distance_norm);
  // event
  fprintf(stream,"%f\t",predictors->first_map_event_distance_norm);
  // sub_event
  fprintf(stream,"%f\t",predictors->subdominant_event_distance_norm);
  // swg
  fprintf(stream,"%f\t",predictors->first_map_swg_score_norm);
  // sub_swg
  fprintf(stream,"%f\t",predictors->subdominant_swg_score_norm);
  // mcs
  fprintf(stream,"%"PRIu64"\t",predictors->mcs);
  // max_region_length
  fprintf(stream,"%f\t",predictors->max_region_length_norm);
  // fs_matches
  fprintf(stream,"%"PRIu64"\t",predictors->first_stratum_matches);
  // sub_matches
  fprintf(stream,"%"PRIu64"\t",predictors->subdominant_stratum_matches);
  // delta
  fprintf(stream,"%"PRIu64"",predictors->subdominant_event_distance-predictors->first_map_event_distance);
}
void matches_predictors_print(
    FILE* const stream,matches_predictors_t* const predictors,
    const char* const read_tag,const uint8_t mapq_score) {
  // Print base predictors
  matches_predictors_print_basic_fields(stream,predictors,read_tag,mapq_score);
  // sub_cand
  fprintf(stream,"\t%"PRIu64"\n",predictors->subdominant_candidates_end1);
}
void paired_matches_predictors_print(
    FILE* const stream,matches_predictors_t* const predictors,
    const char* const read_tag,const uint8_t mapq_score_pair,
    const uint8_t mapq_score_end1,const uint8_t mapq_score_end2) {
  // Print base predictors
  matches_predictors_print_basic_fields(stream,predictors,read_tag,mapq_score_pair);
  // sub_cand_end1
  fprintf(stream,"\t%"PRIu64"\t",predictors->subdominant_candidates_end1);
  // sub_cand_end2
  fprintf(stream,"%"PRIu64"\t",predictors->subdominant_candidates_end2);
  // sigma
  fprintf(stream,"%f\t",predictors->first_map_template_size_sigma);
  // sub_sigma
  fprintf(stream,"%f\t",predictors->subdominant_template_size_sigma);
  // mapq_score_end1
  fprintf(stream,"%d\t",mapq_score_end1);
  // mapq_score_end2
  fprintf(stream,"%d\n",mapq_score_end2);
}

