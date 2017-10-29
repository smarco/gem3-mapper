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
#include "matches/classify/matches_classify.h"
#include "matches/classify/paired_matches_classify.h"
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
    const double primary_error_quality,
    const double primary_template_size_sigma) {
  // Parameters
  const uint64_t read_length = metrics->read_length;
  const int32_t swg_match_score = metrics->swg_match_score;
  const double swg_norm_factor = swg_match_score*read_length;
  // Primary Match
  predictors->primary_edit_distance_norm = DISTANCE_NORMALIZE(primary_edit_distance,read_length);
  predictors->primary_event_distance_norm = DISTANCE_NORMALIZE(primary_event_distance,read_length);
  predictors->primary_swg_score_norm = SCORE_NORMALIZE(primary_swg_score,swg_norm_factor);
  predictors->primary_error_quality = primary_error_quality;
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
    const double primary_error_quality,
    const double primary_template_size_sigma) {
  // Parameters
  const uint64_t read_length = metrics->read_length;
  const int32_t swg_match_score = metrics->swg_match_score;
  const double swg_norm_factor = swg_match_score*read_length;
  // Primary Match
  predictors->primary_edit_distance_norm = DISTANCE_NORMALIZE(primary_edit_distance,read_length);
  predictors->primary_event_distance_norm = DISTANCE_NORMALIZE(primary_event_distance,read_length);
  predictors->primary_swg_score_norm = SCORE_NORMALIZE(primary_swg_score,swg_norm_factor);
  predictors->primary_error_quality = primary_error_quality;
  predictors->primary_template_size_sigma_norm = primary_template_size_sigma/MAX_TEMPLATE_LENGTH_SIGMAS;
  // Subdominant Match
  predictors->subdominant_edit_distance_norm = DISTANCE_NORMALIZE(metrics->min2_edit_distance,read_length);
  predictors->subdominant_event_distance_norm = DISTANCE_NORMALIZE(metrics->min2_event_distance,read_length);
  predictors->subdominant_swg_score_norm = SCORE_NORMALIZE(metrics->max2_swg_score,swg_norm_factor);
  predictors->subdominant_template_size_sigma_norm = 0.0;
}
/*
 * Matches Compute Predictors
 */
void matches_predictors_compute_se(
    matches_predictors_t* const predictors,
    matches_t* const matches) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches==0) return;
  // Classify
  matches_classify(matches);
  // Search Scope
  predictors->candidates_subdominant = matches_counters_count_first_subdominant(matches->counters);
  predictors->candidates_accepted = metrics->candidates_accepted;
  predictors->mcs_end1 = matches->max_complete_stratum;
  predictors->mcs_end2 = 0;
  // Primary/Subdominant predictors
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  match_trace_t* const primary_match = match_traces[0];
  if (num_matches==1) {
    matches_predictors_compute_unique(
        predictors,metrics,primary_match->edit_distance,
        primary_match->event_distance,primary_match->swg_score,
        primary_match->error_quality,MAX_TEMPLATE_LENGTH_SIGMAS);
  } else {
    matches_predictors_compute_mmaps(
        predictors,metrics,primary_match->edit_distance,
        primary_match->event_distance,primary_match->swg_score,
        primary_match->error_quality,MAX_TEMPLATE_LENGTH_SIGMAS);
  }
  // Mappability
  predictors->sequence_length_norm = (double)metrics->read_length/metrics->proper_length;
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
  const uint64_t num_pairs = paired_matches_get_num_maps(paired_matches);
  if (num_pairs==0) return;
  // Classify
  paired_matches_classify(paired_matches);
  // Search Scope
  predictors->candidates_subdominant = matches_counters_count_first_subdominant(paired_matches->counters);
  predictors->candidates_accepted = metrics->candidates_accepted;
  predictors->mcs_end1 = paired_matches_get_max_complete_stratum(paired_matches);
  predictors->mcs_end2 = 0;
  // Primary/Subdominant predictors
  paired_map_t* const primary_match = paired_matches_get_primary_map(paired_matches);
  if (num_pairs==1) {
    matches_predictors_compute_unique(
        predictors,metrics,primary_match->edit_distance,
        primary_match->event_distance,primary_match->swg_score,
        primary_match->error_quality,primary_match->template_length_sigma);
  } else {
    matches_predictors_compute_mmaps(
        predictors,metrics,primary_match->edit_distance,
        primary_match->event_distance,primary_match->swg_score,
        primary_match->error_quality,primary_match->template_length_sigma);
  }
  // Mappability
  predictors->sequence_length_norm = (double)metrics->read_length/metrics->proper_length;
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
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification) {
  /*
   * HEADER
   */
  fprintf(stream,
      "%s" MP_SEP             /* tag   */
      "%02"PRId64 MP_SEP      /* delta */
      "%02"PRId64 MP_SEP      /* wdelta */
      MP_DOUBLE_FORMAT MP_SEP /* edit   */
      MP_DOUBLE_FORMAT MP_SEP /* event  */
      MP_DOUBLE_FORMAT MP_SEP /* swg    */
      MP_DOUBLE_FORMAT MP_SEP /* equal  */
      MP_DOUBLE_FORMAT MP_SEP /* sigma  */
      MP_DOUBLE_FORMAT MP_SEP /* mapq_1 */
      MP_DOUBLE_FORMAT MP_SEP /* mapq_2 */
      MP_DOUBLE_FORMAT MP_SEP /* sub_edit  */
      MP_DOUBLE_FORMAT MP_SEP /* sub_event */
      MP_DOUBLE_FORMAT MP_SEP /* sub_swg   */
      MP_DOUBLE_FORMAT MP_SEP /* sub_sigma */
      "%02"PRIu64 MP_SEP      /* scand  */
      "%02"PRIu64 MP_SEP      /* acand  */
      "%02"PRIu64 MP_SEP      /* mcs1   */
      "%02"PRIu64 MP_SEP      /* mcs2   */
      MP_DOUBLE_FORMAT MP_SEP /* len    */
      MP_DOUBLE_FORMAT MP_SEP /* arlen  */
      MP_DOUBLE_FORMAT MP_SEP /* mrlen  */
      MP_DOUBLE_FORMAT "\n",  /* kmer   */
      /*
       * tag class delta wdelta
       */
      sequence_tag,                          /* tag   */
      matches_classification->delta_group,   /* delta */
      matches_classification->wdelta_group,  /* wdelta */
      /*
       * edit event swg equal sigma mapq_1 mapq_2
       */
      matches_predictors->primary_edit_distance_norm,       /* edit   */
      matches_predictors->primary_event_distance_norm,      /* event  */
      matches_predictors->primary_swg_score_norm,           /* swg    */
      matches_predictors->primary_error_quality,            /* equal  */
      matches_predictors->primary_template_size_sigma_norm, /* sigma  */
      (double)matches_predictors->mapq_end1/60.0,           /* mapq_1 */
      (double)matches_predictors->mapq_end2/60.0,           /* mapq_2 */
      /*
       * sub_edit sub_event sub_swg sub_equal sub_sigma scand acand
       */
      matches_predictors->subdominant_edit_distance_norm,       /* sub_edit  */
      matches_predictors->subdominant_event_distance_norm,      /* sub_event */
      matches_predictors->subdominant_swg_score_norm,           /* sub_swg   */
      matches_predictors->subdominant_template_size_sigma_norm, /* sub_sigma */
      matches_predictors->candidates_subdominant,               /* scand     */
      matches_predictors->candidates_accepted,                  /* acand     */
      /*
       * mcs1 mcs2 len arlen mrlen kmer
       */
      matches_predictors->mcs_end1,               /* mcs1   */
      matches_predictors->mcs_end2,               /* mcs2   */
      matches_predictors->sequence_length_norm,   /* len    */
      matches_predictors->avg_region_length_norm, /* arlen  */
      matches_predictors->max_region_length_norm, /* mrlen  */
      matches_predictors->kmer_frequency);        /* kmer   */
}
void matches_predictors_se_print(
    FILE* const stream,
    const char* const sequence_tag,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification) {
  // Unmapped
  if (matches_classification->matches_class==matches_class_unmapped) return;
  // Print predictors
  matches_predictors_print(stream,sequence_tag,matches_predictors,matches_classification);
}
void matches_predictors_pe_print(
    FILE* const stream,
    const char* const sequence_tag,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification) {
  // Unmapped
  if (matches_classification->matches_class==matches_class_unmapped) return;
  // Print predictors
  matches_predictors_print(stream,sequence_tag,matches_predictors,matches_classification);
}
