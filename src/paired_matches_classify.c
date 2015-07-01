/*
 * PROJECT: GEMMapper
 * FILE: paired_matches_classify.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "paired_matches_classify.h"

/*
 * PE Classify
 */
GEM_INLINE matches_class_t paired_matches_classify(paired_matches_t* const paired_matches) {
  // Parameters
  matches_metrics_t* const metrics = &paired_matches->metrics;
  // Classify
  if (metrics->total_matches_sampled == 0) {
    return matches_class_unmapped;
  } else if (metrics->total_matches_sampled == 1) {
    return matches_class_unique;
  } else {
    matches_t* const matches_end1 = paired_matches->matches_end1;
    matches_t* const matches_end2 = paired_matches->matches_end2;
    paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
    match_trace_t* const primary_end1 = paired_map_get_match_end1(paired_matches,paired_map);
    match_trace_t* const primary_end2 = paired_map_get_match_end2(paired_matches,paired_map);
    match_trace_t* const subdominant_end1 = paired_map_get_match_end1(paired_matches,paired_map+1);
    match_trace_t* const subdominant_end2 = paired_map_get_match_end2(paired_matches,paired_map+1);
    if (matches_cigar_cmp(matches_end1->cigar_vector,primary_end1,matches_end1->cigar_vector,subdominant_end1)==0 &&
        matches_cigar_cmp(matches_end2->cigar_vector,primary_end2,matches_end2->cigar_vector,subdominant_end2)==0) {
      return matches_class_tie_indistinguishable;
    } else if (metrics->max1_swg_score == metrics->max2_swg_score) {
      return matches_class_tie_swg_score;
    } else if (metrics->min1_edit_distance == metrics->min2_edit_distance) {
      return matches_class_tie_edit_distance;
    } else if (metrics->min1_counter_value == metrics->min2_counter_value ||
               metrics->min1_counter_value+1 == metrics->min2_counter_value) {
      return matches_class_tie_event_distance;
    } else {
      return matches_class_mmap;
    }
  }
}
GEM_INLINE void paired_matches_classify_compute_predictors(
    paired_matches_t* const paired_matches,matches_predictors_t* const predictors,
    const swg_penalties_t* const swg_penalties,const uint64_t total_read_length,
    const uint64_t max_region_length,uint64_t const proper_length,
    uint64_t const overriding_mcs) {
  // Compute SE predictors
  const uint64_t mcs = (overriding_mcs!=UINT64_MAX) ? overriding_mcs : paired_matches->max_complete_stratum;
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches==0) {
    matches_classify_compute_predictors_unmapped(predictors,&paired_matches->metrics,
        total_read_length,max_region_length,proper_length,mcs);
    // Compute PE specific predictors
    predictors->first_map_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
    predictors->subdominant_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
  } else {
    paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
    const uint64_t primary_map_distance = paired_map[0].distance;
    const uint64_t primary_map_edit_distance = paired_map[0].edit_distance;
    const int32_t primary_map_swg_score = paired_map[0].swg_score;
    matches_classify_compute_predictors_mapped(predictors,&paired_matches->metrics,
        primary_map_distance,primary_map_edit_distance,primary_map_swg_score,
        swg_penalties,total_read_length,max_region_length,proper_length,mcs);
    // Compute PE specific predictors
    predictors->first_map_template_size_sigma = paired_matches->metrics.min1_template_length_sigma;
    predictors->subdominant_template_size_sigma = paired_matches->metrics.min2_template_length_sigma;
    match_trace_t* const match_end1 = paired_map_get_match_end1(paired_matches,paired_map);
    match_trace_t* const match_end2 = paired_map_get_match_end1(paired_matches,paired_map);
    predictors->mapq_end1 = match_end1->mapq_score;
    predictors->mapq_end2 = match_end2->mapq_score;
  }
  // First/Subdominant strata
  predictors->first_stratum_matches = paired_matches_get_first_stratum_matches(paired_matches);
  predictors->subdominant_stratum_matches = paired_matches_get_subdominant_stratum_matches(paired_matches);
  // Subdominant candidates
  predictors->subdominant_candidates_end1 = paired_matches->matches_end1->metrics.subdominant_candidates;
  predictors->subdominant_candidates_end2 = paired_matches->matches_end2->metrics.subdominant_candidates;
}
GEM_INLINE double paired_matches_classify_unique(matches_predictors_t* const predictors) {
  const double lr_factor = -343.3571 +
      (double)predictors->first_map_edit_distance_norm * 46.6914 +
      (double)predictors->first_map_event_distance_norm * 301.8511 +
      (double)predictors->mcs * 2.3967 +
      (double)predictors->max_region_length_norm * -0.5880 +
      (double)predictors->first_map_template_size_sigma * -0.5061 +
      (double)predictors->mapq_end1 * 0.1308;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
GEM_INLINE double paired_matches_classify_mmaps(matches_predictors_t* const predictors) {
  const double lr_factor = -421.36852 +
      (double)predictors->first_map_edit_distance_norm * 75.46711 +
      (double)predictors->subdominant_edit_distance_norm * -54.84610 +
      (double)predictors->first_map_event_distance_norm * 382.99138 +
      (double)predictors->first_map_swg_score_norm * 18.25120 +
      (double)predictors->mcs * 3.13206 +
      (double)predictors->max_region_length_norm * -0.21534 +
      (double)predictors->mapq_end1 * 0.06624;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
GEM_INLINE double paired_matches_classify_ties(matches_predictors_t* const predictors) {
  const double lr_factor = -2.093e+02 +
      (double)predictors->first_map_edit_distance_norm * 1.465e+02 +
      (double)predictors->subdominant_edit_distance_norm * -1.316e+02 +
      (double)predictors->first_map_event_distance_norm * 6.239e+02 +
      (double)predictors->subdominant_event_distance_norm * -4.275e+02 +
      (double)predictors->first_map_swg_score_norm * 4.753e+01 +
      (double)predictors->subdominant_swg_score_norm * -5.058e+01 +
      (double)predictors->mcs * 5.029e-01 +
      (double)predictors->max_region_length_norm * -1.069e-01 +
      (double)predictors->subdominant_stratum_matches * -1.142e-03 +
      (double)predictors->first_map_template_size_sigma * -8.142e-01 +
      (double)predictors->subdominant_template_size_sigma * 7.452e-01 +
      (double)predictors->mapq_end1 * 1.818e-02;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
