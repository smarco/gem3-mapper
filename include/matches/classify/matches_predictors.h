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

#ifndef MATCHES_PREDICTORS_H_
#define MATCHES_PREDICTORS_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/paired_matches.h"
#include "matches/classify/matches_classify.h"

/*
 * Predictors
 */
typedef struct {
  /* Primary Match */
  double primary_edit_distance_norm;
  double primary_event_distance_norm;
  double primary_swg_score_norm;
  double primary_error_quality;
  double primary_template_size_sigma_norm;
  /* Subdominant Match */
  double subdominant_edit_distance_norm;
  double subdominant_event_distance_norm;
  double subdominant_swg_score_norm;
  double subdominant_template_size_sigma_norm;
  uint64_t candidates_subdominant;
  uint64_t candidates_accepted;
  /* Search Scope */
  uint64_t mcs_end1;
  uint64_t mcs_end2;
  /* Mappability */
  double sequence_length_norm;
  double avg_region_length_norm;
  double max_region_length_norm;
  double kmer_frequency;
  /* MAPQ */
  uint8_t mapq_end1;
  uint8_t mapq_end2;
} matches_predictors_t;

/*
 * Matches Compute Predictors
 */
void matches_predictors_compute_se(
    matches_predictors_t* const predictors,
    matches_t* const matches);
void matches_predictors_compute_pe(
    matches_predictors_t* const predictors,
    paired_matches_t* const paired_matches);

/*
 * Display
 */
void matches_predictors_se_print(
    FILE* const stream,
    const char* const sequence_tag,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification);
void matches_predictors_pe_print(
    FILE* const stream,
    const char* const sequence_tag,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification);

#endif /* MATCHES_PREDICTORS_H_ */
