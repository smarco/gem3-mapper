/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "matches/classify/matches_predictors.h"
#include "archive/score/archive_score_se.h"

/*
 * Utils
 */
#define DISTANCE_NORMALIZE(distance,read_length) ((double)(read_length-distance))/((double)read_length)
#define SCORE_NORMALIZE(score,max_score) ((double)score/(double)max_score)

/*
 * Compute Predictors
 */
void matches_predictors_compute_unique(
    matches_predictors_t* const predictors,
    matches_metrics_t* const metrics,
    const uint64_t primary_edit_distance,
    const uint64_t primary_event_distance,
    const int32_t primary_swg_score,
    const double primary_template_size_sigma) {
  // Parameters
  const uint64_t read_length = metrics->read_length;
  const int32_t swg_match_score = metrics->swg_match_score;
  const double swg_norm_factor = swg_match_score*read_length;
  // Primary Match
  predictors->primary_edit_distance_norm = DISTANCE_NORMALIZE(primary_edit_distance,read_length);
  predictors->primary_event_distance_norm = DISTANCE_NORMALIZE(primary_event_distance,read_length);
  predictors->primary_swg_score_norm = SCORE_NORMALIZE(primary_swg_score,swg_norm_factor);
  predictors->primary_template_size_sigma_norm = primary_template_size_sigma/MAX_TEMPLATE_LENGTH_SIGMAS;
  // Subdominant Match
  predictors->subdominant_edit_distance_norm = 0.0;
  predictors->subdominant_event_distance_norm = 0.0;
  predictors->subdominant_swg_score_norm = 0.0;
  predictors->subdominant_template_size_sigma_norm = 1.0;
}
void matches_predictors_compute_mmaps(
    matches_predictors_t* const predictors,
    matches_metrics_t* const metrics,
    const uint64_t primary_edit_distance,
    const uint64_t primary_event_distance,
    const int32_t primary_swg_score,
    const double primary_template_size_sigma,
    const uint64_t subdominant_edit_distance,
    const uint64_t subdominant_event_distance,
    const int32_t subdominant_swg_score,
    const double subdominant_template_size_sigma) {
  // Parameters
  const uint64_t read_length = metrics->read_length;
  const int32_t swg_match_score = metrics->swg_match_score;
  const double swg_norm_factor = swg_match_score*read_length;
  // Primary Match
  predictors->primary_edit_distance_norm = DISTANCE_NORMALIZE(primary_edit_distance,read_length);
  predictors->primary_event_distance_norm = DISTANCE_NORMALIZE(primary_event_distance,read_length);
  predictors->primary_swg_score_norm = SCORE_NORMALIZE(primary_swg_score,swg_norm_factor);
  predictors->primary_template_size_sigma_norm = primary_template_size_sigma/MAX_TEMPLATE_LENGTH_SIGMAS;
  // Subdominant Match
  if (metrics->min_edit_distance_count > 1 ||
      primary_edit_distance > metrics->min_edit_distance) {
    predictors->subdominant_edit_distance_norm = DISTANCE_NORMALIZE(metrics->min_edit_distance,read_length);
  } else {
    predictors->subdominant_edit_distance_norm = DISTANCE_NORMALIZE(subdominant_edit_distance,read_length);
  }
  // Event distance
  if (metrics->min_event_distance_count > 1 ||
      primary_event_distance > metrics->min_event_distance) {
    predictors->subdominant_event_distance_norm = DISTANCE_NORMALIZE(metrics->min_event_distance,read_length);
  } else {
    predictors->subdominant_event_distance_norm = DISTANCE_NORMALIZE(subdominant_event_distance,read_length);
  }
  // SWG
  if (metrics->max_swg_score_count > 1 ||
      primary_swg_score < metrics->max_swg_score) {
    predictors->subdominant_swg_score_norm = SCORE_NORMALIZE(metrics->max_swg_score,swg_norm_factor);
  } else {
    predictors->subdominant_swg_score_norm = SCORE_NORMALIZE(subdominant_swg_score,swg_norm_factor);
  }
  // Template length
  if (metrics->min_template_length_sigma_count > 1 ||
      primary_template_size_sigma > metrics->min_template_length_sigma) {
    predictors->subdominant_template_size_sigma_norm = metrics->min_template_length_sigma/MAX_TEMPLATE_LENGTH_SIGMAS;
  } else {
    predictors->subdominant_template_size_sigma_norm = subdominant_template_size_sigma/MAX_TEMPLATE_LENGTH_SIGMAS;
  }
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
  if (metrics->accepted_matches==0) return;
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  match_trace_t* const primary_match = match_traces[0];
  if (metrics->accepted_matches==1) {
    matches_predictors_compute_unique(
        predictors,metrics,primary_match->edit_distance,
        primary_match->event_distance,primary_match->swg_score,
        MAX_TEMPLATE_LENGTH_SIGMAS);
  } else {
    match_trace_t* const subdominant_match = match_traces[1];
    matches_predictors_compute_mmaps(predictors,metrics,
        primary_match->edit_distance,primary_match->event_distance,
        primary_match->swg_score,MAX_TEMPLATE_LENGTH_SIGMAS,
        subdominant_match->edit_distance,subdominant_match->event_distance,
        subdominant_match->swg_score,MAX_TEMPLATE_LENGTH_SIGMAS);
  }
  // Search Scope
  const double log4 = (double)gem_loge(4);
  predictors->accepted_candidates = (double)gem_loge((float)metrics->accepted_candidates)/log4;
  predictors->aligned_candidates = (double)gem_loge((float)metrics->aligned_alignments)/log4;
  predictors->accepted_matches = (double)gem_loge((float)metrics->accepted_matches)/log4;
  predictors->mcs_end1 = (matches->max_complete_stratum!=ALL) ? matches->max_complete_stratum : 0;
  predictors->mcs_end2 = 0;
  // Mappability
  predictors->avg_region_length_norm = (double)metrics->avg_region_length/metrics->proper_length;
  predictors->max_region_length_norm = (double)metrics->max_region_length/metrics->proper_length;
  predictors->kmer_frequency = metrics->kmer_frequency;
  // MAPQ
  predictors->mapq_end1 = 0;
  predictors->mapq_end2 = 0;
}
/*
 * PE-Matches Compute Predictors
 */
void matches_predictors_compute_pe(
    matches_predictors_t* const predictors,
    paired_matches_t* const paired_matches) {
  // Parameters
  matches_metrics_t* const metrics = &paired_matches->metrics;
  // Primary/Subdominant predictors
  if (metrics->accepted_matches==0) return;
  paired_map_t* const primary_match = paired_matches_get_primary_map(paired_matches);
  if (metrics->accepted_matches==1) {
    matches_predictors_compute_unique(
        predictors,metrics,primary_match->edit_distance,
        primary_match->event_distance,primary_match->swg_score,
        primary_match->template_length_sigma);
  } else {
    paired_map_t* const subdominant_match = paired_matches_get_secondary_map(paired_matches);
    matches_predictors_compute_mmaps(predictors,metrics,
        primary_match->edit_distance,primary_match->event_distance,
        primary_match->swg_score,primary_match->template_length_sigma,
        subdominant_match->edit_distance,subdominant_match->event_distance,
        subdominant_match->swg_score,subdominant_match->template_length_sigma);
  }
  // Search Scope
  const double log4 = (double)gem_loge(4);
  predictors->accepted_candidates = (double)gem_loge((float)metrics->accepted_candidates)/log4;
  predictors->aligned_candidates = (double)gem_loge((float)metrics->aligned_alignments)/log4;
  predictors->accepted_matches = (double)gem_loge((float)metrics->accepted_matches)/log4;
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  predictors->mcs_end1 = (matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0;
  predictors->mcs_end2 = (matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0;
  // Mappability
  predictors->avg_region_length_norm = (double)metrics->avg_region_length/metrics->proper_length;
  predictors->max_region_length_norm = (double)metrics->max_region_length/metrics->proper_length;
  predictors->kmer_frequency = metrics->kmer_frequency;
  // MAPQ
  match_trace_t* const primary_end1 = primary_match->match_trace_end1;
  match_trace_t* const primary_end2 = primary_match->match_trace_end2;
  predictors->mapq_end1 = primary_end1->mapq_score;
  predictors->mapq_end2 = primary_end2->mapq_score;
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
    matches_predictors_t* const matches_predictors) {
  /*
   * HEADER
   *   tp tag class
   *   edit event swg sigma mapq_1 mapq_2
   *   sub_edit sub_event sub_swg sub_sigma
   *   ccand cmatch mcs1 mcs2 mrl kmerf
   */
  fprintf(stream,
      /*
       * tag class
       */
      "%s" MP_SEP /* tag   */
      "%s" MP_SEP /* class */
      /*
       * edit event swg sigma mapq_1 mapq_2
       */
      MP_DOUBLE_FORMAT MP_SEP /* edit   */
      MP_DOUBLE_FORMAT MP_SEP /* event  */
      MP_DOUBLE_FORMAT MP_SEP /* swg    */
      MP_DOUBLE_FORMAT MP_SEP /* sigma  */
      MP_DOUBLE_FORMAT MP_SEP /* mapq_1 */
      MP_DOUBLE_FORMAT MP_SEP /* mapq_2 */
      /*
       * sub_edit sub_event sub_swg sub_sigma
       */
      MP_DOUBLE_FORMAT MP_SEP /* sub_edit  */
      MP_DOUBLE_FORMAT MP_SEP /* sub_event */
      MP_DOUBLE_FORMAT MP_SEP /* sub_swg   */
      MP_DOUBLE_FORMAT MP_SEP /* sub_sigma */
      /*
       * ccand acand cmatch mcs1 mcs2 arl mrl kmerf
       */
      MP_DOUBLE_FORMAT MP_SEP /* ccand  */
      MP_DOUBLE_FORMAT MP_SEP /* acand  */
      MP_DOUBLE_FORMAT MP_SEP /* cmatch */
      "%02"PRIu64 MP_SEP      /* mcs1   */
      "%02"PRIu64 MP_SEP      /* mcs2   */
      MP_DOUBLE_FORMAT MP_SEP /* arl    */
      MP_DOUBLE_FORMAT MP_SEP /* mrl    */
      MP_DOUBLE_FORMAT "\n",  /* kmerf  */
      /*
       * tag class
       */
      sequence_tag, /* tag   */
      match_class,  /* class */
      /*
       * edit event swg sigma mapq_1 mapq_2
       */
      matches_predictors->primary_edit_distance_norm,       /* edit   */
      matches_predictors->primary_event_distance_norm,      /* event  */
      matches_predictors->primary_swg_score_norm,           /* swg    */
      matches_predictors->primary_template_size_sigma_norm, /* sigma  */
      (double)matches_predictors->mapq_end1/60.0,           /* mapq_1 */
      (double)matches_predictors->mapq_end2/60.0,           /* mapq_2 */
      /*
       * sub_edit sub_event sub_swg sub_sigma
       */
      matches_predictors->subdominant_edit_distance_norm,       /* sub_edit  */
      matches_predictors->subdominant_event_distance_norm,      /* sub_event */
      matches_predictors->subdominant_swg_score_norm,           /* sub_swg   */
      matches_predictors->subdominant_template_size_sigma_norm, /* sub_sigma */
      /*
       * ccand acand cmatch mcs1 mcs2 arl mrl kmerf
       */
      matches_predictors->accepted_candidates,    /* ccand  */
      matches_predictors->aligned_candidates,     /* acand  */
      matches_predictors->accepted_matches,       /* cmatch */
      matches_predictors->mcs_end1,               /* mcs1   */
      matches_predictors->mcs_end2,               /* mcs2   */
      matches_predictors->avg_region_length_norm, /* arl    */
      matches_predictors->max_region_length_norm, /* mrl    */
      matches_predictors->kmer_frequency);        /* kmerf  */
}
void matches_predictors_se_print(
    FILE* const stream,
    const char* const sequence_tag,
    const matches_class_t matches_class,
    matches_predictors_t* const matches_predictors) {
  // Class
  switch (matches_class) {
    case matches_class_unmapped:
      break;
    case matches_class_tie_perfect:
    case matches_class_tie:
      matches_predictors_print(stream,sequence_tag,"tie   ",matches_predictors);
      break;
    case matches_class_mmap_d1:
      matches_predictors_print(stream,sequence_tag,"mmapD1",matches_predictors);
      break;
    case matches_class_mmap:
      matches_predictors_print(stream,sequence_tag,"mmap  ",matches_predictors);
      break;
    case matches_class_unique:
      matches_predictors_print(stream,sequence_tag,"unique",matches_predictors);
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
    matches_predictors_t* const matches_predictors) {
  // Unmapped
  if (paired_matches_class==paired_matches_class_unmapped) return;
  // High quality ends
  const bool high_quality_ends =
      matches_predictors->mapq_end1 >= MAPQ_CONFIDENCE_SCORE_MIN &&
      matches_predictors->mapq_end2 >= MAPQ_CONFIDENCE_SCORE_MIN;
  if (high_quality_ends) {
    matches_predictors_print(stream,sequence_tag,"signal ",matches_predictors);
    return;
  }
  switch (paired_matches_class) {
    case paired_matches_class_unmapped:
      break;
    case paired_matches_class_tie_perfect:
    case paired_matches_class_tie:
      matches_predictors_print(stream,sequence_tag,"tie ",matches_predictors);
      break;
    case paired_matches_class_mmap_d1:
      matches_predictors_print(stream,sequence_tag,"mmapD1",matches_predictors);
      break;
    case paired_matches_class_mmap:
      matches_predictors_print(stream,sequence_tag,"mmap  ",matches_predictors);
      break;
    case paired_matches_class_unique:
      matches_predictors_print(stream,sequence_tag,"unique",matches_predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
