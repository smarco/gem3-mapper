/*
 * PROJECT: GEMMapper
 * FILE: matches_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches/matches_predictors.h"

/*
 * Compute Predictors
 */
void matches_predictors_compute_unmapped(matches_predictors_t* const predictors,const uint64_t read_length) {
  // Compute event distances
  predictors->first_map_event_distance = read_length;
  predictors->subdominant_event_distance = read_length;
  predictors->first_map_edit_distance = read_length;
  predictors->subdominant_edit_distance = read_length;
  predictors->first_map_swg_score = 0;
  predictors->subdominant_swg_score = 0;
  // Normalize event distances
  predictors->first_map_event_distance_norm = 0.0;
  predictors->subdominant_event_distance_norm = 0.0;
  // Normalize edit distances
  predictors->first_map_edit_distance_norm = 0.0;
  predictors->subdominant_edit_distance_norm = 0.0;
  // Normalize SWG scores
  predictors->first_map_swg_score_norm =  0.0;
  predictors->subdominant_swg_score_norm =  0.0;
}
void matches_predictors_compute_mapped(
    matches_predictors_t* const predictors,
    matches_metrics_t* const matches_metrics,
    const uint64_t primary_map_distance,
    const uint64_t primary_map_edit_distance,
    const int32_t primary_map_swg_score,
    const uint64_t read_length,
    const int32_t swg_match_score) {
  // Compute primary/subdominant event distance
  predictors->first_map_event_distance = primary_map_distance;
  if (primary_map_distance!=matches_metrics->min1_counter_value) {
    predictors->subdominant_event_distance = matches_metrics->min1_counter_value;
  } else {
    predictors->subdominant_event_distance = matches_metrics->min2_counter_value;
  }
  if (predictors->subdominant_event_distance==UINT32_MAX) {
    predictors->subdominant_event_distance = read_length;
  }
  // Compute primary/subdominant edit distance
  predictors->first_map_edit_distance = primary_map_edit_distance;
  if (primary_map_edit_distance!=matches_metrics->min1_edit_distance) {
    predictors->subdominant_edit_distance = matches_metrics->min1_edit_distance;
  } else {
    predictors->subdominant_edit_distance = matches_metrics->min2_edit_distance;
  }
  if (predictors->subdominant_edit_distance==UINT32_MAX) {
    predictors->subdominant_edit_distance = read_length;
  }
  // Compute primary/subdominant SWG score
  predictors->first_map_swg_score = primary_map_swg_score;
  if (predictors->first_map_swg_score < 0) predictors->first_map_swg_score = 0;
  if (primary_map_swg_score!=matches_metrics->max1_swg_score) {
    predictors->subdominant_swg_score = matches_metrics->max1_swg_score;
  } else {
    predictors->subdominant_swg_score = matches_metrics->max2_swg_score;
  }
  if (predictors->subdominant_swg_score < 0) predictors->subdominant_swg_score = 0;
  // Normalize event distances
  predictors->first_map_event_distance_norm =
      ((double)(read_length-predictors->first_map_event_distance))/((double)read_length);
  predictors->subdominant_event_distance_norm =
      ((double)(read_length-predictors->subdominant_event_distance))/((double)read_length);
  // Normalize edit distances
  predictors->first_map_edit_distance_norm =
      ((double)(read_length-predictors->first_map_edit_distance))/((double)read_length);
  predictors->subdominant_edit_distance_norm =
      ((double)(read_length-predictors->subdominant_edit_distance))/((double)read_length);
  // Normalize SWG scores
  const double swg_norm_factor = swg_match_score*read_length;
  predictors->first_map_swg_score_norm = (double)predictors->first_map_swg_score/swg_norm_factor;
  predictors->subdominant_swg_score_norm = (double)predictors->subdominant_swg_score/swg_norm_factor;
}
/*
 * SE Compute Predictors
 */
void matches_predictors_compute(
    matches_t* const matches,
    matches_predictors_t* const predictors,
    approximate_search_metrics_t* const search_metrics,
    const uint64_t mcs) {
  // Parameters
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  matches_metrics_t* const metrics = &matches->metrics;
  // Init
  if (num_matches==0) {
    // First stratum
    predictors->first_stratum_matches = 0;
    // Primary/Subdominant predictors
    matches_predictors_compute_unmapped(predictors,search_metrics->read_length);
  } else {
    // First stratum
    predictors->first_stratum_matches = matches_get_first_stratum_matches(matches);
    // Primary/Subdominant predictors
    match_trace_t* const match = matches_get_match_trace_buffer(matches);
    matches_predictors_compute_mapped(predictors,metrics,
        match->distance,match->edit_distance,match->swg_score,
        search_metrics->read_length,search_metrics->swg_match_score);
  }
  /*
   * Search Scope
   */
  predictors->mcs_end1 = mcs;
  predictors->mcs_end2 = 0;
  predictors->accepted_candidates_end1 = matches->metrics.accepted_candidates;
  predictors->accepted_candidates_end2 = 0;
  /*
   * Mappability
   */
  predictors->max_region_length_norm = (double)search_metrics->max_region_length/search_metrics->proper_length;
  predictors->mappability_p = search_metrics->mappability_p;
  predictors->mappability_2p = search_metrics->mappability_2p;
  /*
   * Template Size
   */
  predictors->first_map_template_size_sigma = 0;
  predictors->subdominant_template_size_sigma = 0;
  /*
   * MAPQ Score
   */
  predictors->mapq_end1 = 0;
  predictors->mapq_end2 = 0;
}
/*
 * PE Compute Predictors
 */
void paired_matches_predictors_compute(
    paired_matches_t* const paired_matches,
    matches_predictors_t* const predictors,
    approximate_search_metrics_t* const search_metrics_end1,
    approximate_search_metrics_t* const search_metrics_end2,
    const uint64_t mcs_end1,
    const uint64_t mcs_end2) {
  // Parameters
  matches_metrics_t* const metrics_end1 = &paired_matches->matches_end1->metrics;
  matches_metrics_t* const metrics_end2 = &paired_matches->matches_end2->metrics;
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  const uint64_t total_read_length = search_metrics_end1->read_length + search_metrics_end2->read_length;
  // Compute PE predictors
  if (num_matches==0) {
    // First stratum
    predictors->first_stratum_matches = 0;
    // Primary/Subdominant predictors
    matches_predictors_compute_unmapped(predictors,total_read_length);
    // Template Size
    predictors->first_map_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
    predictors->subdominant_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
    // MAPQ Score
    predictors->mapq_end1 = 0;
    predictors->mapq_end2 = 0;
  } else {
    // First stratum
    predictors->first_stratum_matches = paired_matches_get_first_stratum_matches(paired_matches);
    paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
    // Primary/Subdominant predictors
    matches_predictors_compute_mapped(predictors,
        &paired_matches->metrics,paired_map[0].distance,
        paired_map[0].edit_distance,paired_map[0].swg_score,
        total_read_length,search_metrics_end1->swg_match_score);
    // Template Size
    predictors->first_map_template_size_sigma = paired_matches->metrics.min1_template_length_sigma;
    predictors->subdominant_template_size_sigma = paired_matches->metrics.min2_template_length_sigma;
    // MAPQ Score
    match_trace_t* const match_end1 = paired_map_get_match_end1(paired_matches,paired_map);
    match_trace_t* const match_end2 = paired_map_get_match_end2(paired_matches,paired_map);
    predictors->mapq_end1 = match_end1->mapq_score;
    predictors->mapq_end2 = match_end2->mapq_score;
  }
  /*
   * Search Scope
   */
  predictors->mcs_end1 = mcs_end1;
  predictors->mcs_end2 = mcs_end2;
  predictors->accepted_candidates_end1 = metrics_end1->accepted_candidates;
  predictors->accepted_candidates_end2 = metrics_end2->accepted_candidates;
  /*
   * Mappability
   */
  const uint64_t max_region_length_end1 = (search_metrics_end1->max_region_length==UINT32_MAX) ?
      0 : search_metrics_end1->max_region_length;
  const uint64_t max_region_length_end2 = (search_metrics_end2->max_region_length==UINT32_MAX) ?
      0 : search_metrics_end2->max_region_length;
  const uint64_t max_region_length = MAX(max_region_length_end1,max_region_length_end2);
  predictors->max_region_length_norm = (double)max_region_length/search_metrics_end1->proper_length;
  predictors->mappability_p = MAX(search_metrics_end1->mappability_p,search_metrics_end2->mappability_p);
  predictors->mappability_2p = MAX(search_metrics_end1->mappability_2p,search_metrics_end2->mappability_2p);
}
#define MP_SEP           "  "
#define MP_DOUBLE_FORMAT "%f"
void matches_predictors_print(
    FILE* const stream,
    const char* const sequence_tag,
    const char* const match_class,
    matches_predictors_t* const predictors) {
  // Tag
  fprintf(stream,"%s" MP_SEP,sequence_tag);
  // Class
  fprintf(stream,"%s" MP_SEP,match_class);
  /*
   * Distance
   */
  // edit
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->first_map_edit_distance_norm);
  // sub_edit
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->subdominant_edit_distance_norm);
  // event
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->first_map_event_distance_norm);
  // sub_event
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->subdominant_event_distance_norm);
  // swg
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->first_map_swg_score_norm);
  // sub_swg
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->subdominant_swg_score_norm);
  /*
   * Search Scope
   */
  // fs_matches
  fprintf(stream,"%03"PRIu64 MP_SEP,predictors->first_stratum_matches);
  // mcs_end1
  fprintf(stream,"%02"PRIu64 MP_SEP,predictors->mcs_end1);
  // mcs_end2
  fprintf(stream,"%02"PRIu64 MP_SEP,predictors->mcs_end2);
  // delta
  fprintf(stream,"%03"PRIu64 MP_SEP,predictors->subdominant_event_distance-predictors->first_map_event_distance);
  // accepted_can_end1
  fprintf(stream,"%03"PRIu64 MP_SEP,predictors->accepted_candidates_end1);
  // accepted_can_end2
  fprintf(stream,"%03"PRIu64 MP_SEP,predictors->accepted_candidates_end2);
  /*
   * Mappability
   */
  // max_region_length
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->max_region_length_norm);
  // Mappability (p)
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->mappability_p);
  // Mappability (2p)
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->mappability_2p);
  /*
   * Template Size
   */
  // sigma
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->first_map_template_size_sigma);
  // sub_sigma
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,predictors->subdominant_template_size_sigma);
  /*
   * MAPQ
   */
  // mapq_1
  fprintf(stream,"%02d" MP_SEP,predictors->mapq_end1);
  // mapq_2
  fprintf(stream,"%02d\n",predictors->mapq_end2);
}
void matches_predictors_se_print(
    FILE* const stream,
    const char* const sequence_tag,
    const matches_class_t matches_class,
    matches_predictors_t* const predictors) {
  // Class
  switch (matches_class) {
    case matches_class_unmapped:
      break;
    case matches_class_tie_d0:
      matches_predictors_print(stream,sequence_tag,"noise ",predictors);
      break;
    case matches_class_tie_d1:
      matches_predictors_print(stream,sequence_tag,"tie   ",predictors);
      break;
    case matches_class_mmap:
      matches_predictors_print(stream,sequence_tag,"mmap  ",predictors);
      break;
    case matches_class_unique:
      matches_predictors_print(stream,sequence_tag,"unique",predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void matches_predictors_pe_print(
    FILE* const stream,
    const char* const sequence_tag,
    const paired_matches_class_t paired_matches_class,
    matches_predictors_t* const predictors) {
  switch (paired_matches_class) {
    case paired_matches_class_unmapped:
      break;
    case paired_matches_class_subdominant_end:
    case paired_matches_class_tie_d0:
      matches_predictors_print(stream,sequence_tag,"noise ",predictors);
      break;
    case paired_matches_class_tie_d1: {
      matches_predictors_print(stream,sequence_tag,"tie   ",predictors);
      break;
    }
    case paired_matches_class_mmap: {
      matches_predictors_print(stream,sequence_tag,"mmap  ",predictors);
      break;
    }
    case paired_matches_class_unique: {
      matches_predictors_print(stream,sequence_tag,"unique",predictors);
      break;
    }
    case paired_matches_class_high_quality_ends:
      matches_predictors_print(stream,sequence_tag,"signal",predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
