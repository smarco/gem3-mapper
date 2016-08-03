/*
 * PROJECT: GEMMapper
 * FILE: matches_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches/matches_predictors.h"

/*
 * Utils
 */
#define DISTANCE_NORMALIZE(distance,read_length) ((double)(read_length-distance))/((double)read_length)
#define SCORE_NORMALIZE(score,max_score) ((double)score/(double)max_score)

/*
 * Match Compute Predictors
 */
void match_predictors_compute_se(
    match_predictors_t* const match_predictors,
    matches_t* const matches,
    match_trace_t* const match) {
  // Unammped
  if (match==NULL) {
    match_predictors->map_edit_distance_norm = 0.0;
    match_predictors->map_event_distance_norm = 0.0;
    match_predictors->map_swg_score_norm = 0.0;
  } else {
    // Parameters
    matches_metrics_t* const metrics = &matches->metrics;
    const uint64_t read_length = metrics->read_length;
    const int32_t swg_match_score = metrics->swg_match_score;
    // Match
    const double swg_norm_factor = (double)swg_match_score*(double)read_length;
    match_predictors->map_edit_distance_norm = DISTANCE_NORMALIZE(match->edit_distance,read_length);
    match_predictors->map_event_distance_norm = DISTANCE_NORMALIZE(match->distance,read_length);
    match_predictors->map_swg_score_norm = SCORE_NORMALIZE(match->swg_score,swg_norm_factor);
  }
  match_predictors->map_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
  match_predictors->mapq_end1 = 0;
  match_predictors->mapq_end2 = 0;
}
void match_predictors_compute_pe(
    match_predictors_t* const match_predictors,
    paired_matches_t* const paired_matches,
    paired_map_t* const paired_map) {
  // Unammped
  if (paired_map==NULL) {
    match_predictors->map_edit_distance_norm = 0.0;
    match_predictors->map_event_distance_norm = 0.0;
    match_predictors->map_swg_score_norm = 0.0;
    match_predictors->map_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
    match_predictors->mapq_end1 = 0;
    match_predictors->mapq_end2 = 0;
  } else {
    // Parameters
    matches_metrics_t* const metrics = &paired_matches->metrics;
    const uint64_t read_length = metrics->read_length;
    const int32_t swg_match_score = metrics->swg_match_score;
    // Match
    const double swg_norm_factor = (double)swg_match_score*(double)read_length;
    match_predictors->map_edit_distance_norm = DISTANCE_NORMALIZE(paired_map->edit_distance,read_length);
    match_predictors->map_event_distance_norm = DISTANCE_NORMALIZE(paired_map->distance,read_length);
    match_predictors->map_swg_score_norm = SCORE_NORMALIZE(paired_map->swg_score,swg_norm_factor);
    // Template Size
    match_predictors->map_template_size_sigma = paired_map->template_length_sigma;
    // MAPQ Score
    match_trace_t* const match_end1 = paired_map_get_match_end1(paired_matches,paired_map);
    match_trace_t* const match_end2 = paired_map_get_match_end2(paired_matches,paired_map);
    match_predictors->mapq_end1 = match_end1->mapq_score;
    match_predictors->mapq_end2 = match_end2->mapq_score;
  }
}
/*
 * Compute Predictors
 */
void matches_predictors_compute_unmapped(
    matches_predictors_t* const predictors,
    matches_metrics_t* const metrics) {
  // Normalize event distances
  predictors->best_map_event_distance_norm = 0.0;
  predictors->subdominant_event_distance_norm = 0.0;
  // Normalize edit distances
  predictors->best_map_edit_distance_norm = 0.0;
  predictors->subdominant_edit_distance_norm = 0.0;
  // Normalize SWG scores
  predictors->best_map_swg_score_norm =  0.0;
  predictors->subdominant_swg_score_norm =  0.0;
  // Template Size
  predictors->best_map_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
  predictors->subdominant_template_size_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
}
void matches_predictors_compute_mapped(
    matches_predictors_t* const predictors,
    matches_metrics_t* const metrics) {
  const uint64_t read_length = metrics->read_length;
  const int32_t swg_match_score = metrics->swg_match_score;
  // Auxiliary
  uint64_t best_map_edit_distance;
  uint64_t best_map_event_distance;
  int32_t best_map_swg_score;
  uint64_t subdominant_edit_distance;
  uint64_t subdominant_event_distance;
  int32_t subdominant_swg_score;
  // Compute best/subdominant edit distance
  best_map_edit_distance = metrics->min1_edit_distance;
  subdominant_edit_distance = metrics->min2_edit_distance;
  if (subdominant_edit_distance==UINT32_MAX) {
    subdominant_edit_distance = read_length;
  }
  // Compute best/subdominant event distance
  best_map_event_distance = metrics->min1_counter_value;
  subdominant_event_distance = metrics->min2_counter_value;
  if (subdominant_event_distance==UINT32_MAX) {
    subdominant_event_distance = read_length;
  }
  // Compute best/subdominant SWG score
  best_map_swg_score = metrics->max1_swg_score;
  if (best_map_swg_score < 0) best_map_swg_score = 0;
  subdominant_swg_score = metrics->max2_swg_score;
  if (subdominant_swg_score < 0) subdominant_swg_score = 0;
  // Normalize event distances
  predictors->best_map_event_distance_norm = DISTANCE_NORMALIZE(best_map_event_distance,read_length);
  predictors->subdominant_event_distance_norm = DISTANCE_NORMALIZE(subdominant_event_distance,read_length);
  // Normalize edit distances
  predictors->best_map_edit_distance_norm = DISTANCE_NORMALIZE(best_map_edit_distance,read_length);
  predictors->subdominant_edit_distance_norm = DISTANCE_NORMALIZE(subdominant_edit_distance,read_length);
  // Normalize SWG scores
  const double swg_norm_factor = swg_match_score*read_length;
  predictors->best_map_swg_score_norm = SCORE_NORMALIZE(best_map_swg_score,swg_norm_factor);
  predictors->subdominant_swg_score_norm = SCORE_NORMALIZE(subdominant_swg_score,swg_norm_factor);
  // Template Size
  predictors->best_map_template_size_sigma = metrics->min1_template_length_sigma;
  predictors->subdominant_template_size_sigma = metrics->min2_template_length_sigma;
}
/*
 * Matches Compute Predictors
 */
void matches_predictors_compute_se(
    matches_predictors_t* const predictors,
    matches_t* const matches) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  // Primary/Subdominant predictors
  if (!matches_is_mapped(matches)) {
    matches_predictors_compute_unmapped(predictors,metrics);
  } else {
    matches_predictors_compute_mapped(predictors,metrics);
  }
  // Search Scope
  predictors->matches_accepted = metrics->total_matches_sampled;
  predictors->mcs_end1 = (matches->max_complete_stratum!=ALL) ? matches->max_complete_stratum : 0;
  predictors->mcs_end2 = 0;
  // Mappability
  predictors->max_region_length_norm = (double)metrics->max_region_length/metrics->proper_length;
  predictors->kmer_frequency = metrics->kmer_frequency;
}
/*
 * PE-Matches Compute Predictors
 */
void matches_predictors_compute_pe(
    matches_predictors_t* const predictors,
    paired_matches_t* const paired_matches) {
  // Parameters
  matches_metrics_t* const metrics = &paired_matches->metrics;
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  // Primary/Subdominant predictors
  if (num_matches==0) {
    matches_predictors_compute_unmapped(predictors,metrics);
  } else {
    matches_predictors_compute_mapped(predictors,metrics);
  }
  // Search Scope
  predictors->matches_accepted = metrics->total_matches_sampled;
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  predictors->mcs_end1 = (matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0;
  predictors->mcs_end2 = (matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0;
  // Mappability
  predictors->max_region_length_norm = (double)metrics->max_region_length/metrics->proper_length;
  predictors->kmer_frequency = metrics->kmer_frequency;
}
/*
 * Display
 */
#define MP_SEP           "  "
#define MP_DOUBLE_FORMAT "%f"
void matches_predictors_print(
    FILE* const stream,
    const char* const sequence_tag,
    const char* const match_class,
    matches_predictors_t* const matches_predictors,
    match_predictors_t* const match_predictors) {
  /*
   * HEADER
   *   tp tag class edit event swg sigma mapq_1 mapq_2
   *   edit1 edit2 event1 event2 swg1 swg2 sigma1 sigma2
   *   ack mcs1 mcs2 mrl kmerf
   */
  // tag
  fprintf(stream,"%s" MP_SEP,sequence_tag);
  // class
  fprintf(stream,"%s" MP_SEP,match_class);
  /*
   * edit event swg sigma mapq_1 mapq_2
   */
  // edit
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,match_predictors->map_edit_distance_norm);
  // event
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,match_predictors->map_event_distance_norm);
  // swg
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,match_predictors->map_swg_score_norm);
  // sigma
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,match_predictors->map_template_size_sigma);
  // mapq_1
  fprintf(stream,"%02d" MP_SEP,match_predictors->mapq_end1);
  // mapq_2
  fprintf(stream,"%02d" MP_SEP,match_predictors->mapq_end2);
  /*
   * edit1 edit2 event1 event2 swg1 swg2
   */
  // edit1
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->best_map_edit_distance_norm);
  // edit2
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->subdominant_edit_distance_norm);
  // event1
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->best_map_event_distance_norm);
  // event2
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->subdominant_event_distance_norm);
  // swg1
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->best_map_swg_score_norm);
  // swg2
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->subdominant_swg_score_norm);
  /*
   * sigma1 sigma2
   */
  // sigma1
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->best_map_template_size_sigma);
  // sigma2
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->subdominant_template_size_sigma);
  /*
   * ack mcs1 mcs2 mrl kmerf
   */
  // ack
  fprintf(stream,"%03"PRIu64 MP_SEP,matches_predictors->matches_accepted);
  // mcs1
  fprintf(stream,"%02"PRIu64 MP_SEP,matches_predictors->mcs_end1);
  // mcs2
  fprintf(stream,"%02"PRIu64 MP_SEP,matches_predictors->mcs_end2);
  // mrl
  fprintf(stream,MP_DOUBLE_FORMAT MP_SEP,matches_predictors->max_region_length_norm);
  // kmerf
  fprintf(stream,MP_DOUBLE_FORMAT "\n",matches_predictors->kmer_frequency);
}
void matches_predictors_se_print(
    FILE* const stream,
    const char* const sequence_tag,
    const matches_class_t matches_class,
    matches_predictors_t* const matches_predictors,
    match_predictors_t* const match_predictors) {
  // Class
  switch (matches_class) {
    case matches_class_unmapped:
      break;
    case matches_class_tie_d0:
      matches_predictors_print(stream,sequence_tag,"noise ",matches_predictors,match_predictors);
      break;
    case matches_class_tie_d1:
      matches_predictors_print(stream,sequence_tag,"tie   ",matches_predictors,match_predictors);
      break;
    case matches_class_mmap:
      matches_predictors_print(stream,sequence_tag,"mmap  ",matches_predictors,match_predictors);
      break;
    case matches_class_unique:
      matches_predictors_print(stream,sequence_tag,"unique",matches_predictors,match_predictors);
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
    matches_predictors_t* const matches_predictors,
    match_predictors_t* const match_predictors) {
  // Unmapped
  if (paired_matches_class==paired_matches_class_unmapped) return;
  // High quality ends
  if (match_predictors->mapq_end1>=30 && match_predictors->mapq_end2>=30) {
    matches_predictors_print(stream,sequence_tag,"signal ",matches_predictors,match_predictors);
    return;
  }
  switch (paired_matches_class) {
    case paired_matches_class_unmapped:
      break;
    case paired_matches_class_tie_d0:
      matches_predictors_print(stream,sequence_tag,"noise ",matches_predictors,match_predictors);
      break;
    case paired_matches_class_tie_d1:
      matches_predictors_print(stream,sequence_tag,"tie   ",matches_predictors,match_predictors);
      break;
    case paired_matches_class_mmap:
      matches_predictors_print(stream,sequence_tag,"mmap  ",matches_predictors,match_predictors);
      break;
    case paired_matches_class_unique:
      matches_predictors_print(stream,sequence_tag,"unique",matches_predictors,match_predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
