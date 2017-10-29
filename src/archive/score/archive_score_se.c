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
 * DESCRIPTION:
 *   Archive MAPQ-scoring module produces MAPQ scores for SE-searches
 */

#include "archive/score/archive_score_se.h"
#include "archive/search/archive_search_se.h"
#include "matches/classify/matches_classify.h"
#include "archive/score/archive_score_logit.h"
#include "archive/score/archive_score_logit_models.h"
#include "profiler/profiler.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Archive Scoring Utils
 */
uint8_t archive_score_probability_scale(
    const double probability,
    const double sum_probability,
    const uint8_t floor,
    const uint8_t ceil) {
  const double range = ceil-floor;
  return floor + (uint8_t)round((range*probability)/sum_probability);
}
uint64_t archive_score_probability_to_mapq(const double probability) {
  return (uint64_t)round(-10.0*log10(1.0-probability));
}
/*
 * Archive Scoring Logit Classes
 */
const archive_score_logit_model_t* archive_score_matches_se_get_model(
    search_parameters_t* const search_parameters) {
  if (search_parameters->mapping_mode==mapping_adaptive_filtering_fast) {
    return &logit_model_single_end_fast;
  } else {
    return &logit_model_single_end_sensitive;
  }
}
uint8_t archive_score_matches_se_logit_unique(
    search_parameters_t* const search_parameters,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification) {
  // Fetch logit model
  const archive_score_logit_model_t* const classify_model =
      archive_score_matches_se_get_model(search_parameters);
  // Compute probability
  const double pr = archive_score_logit_unique(
      classify_model,matches_predictors,matches_classification);
  const uint64_t mapq_score = archive_score_probability_to_mapq(pr);
  return (uint8_t)MIN(mapq_score,59);
}
uint8_t archive_score_matches_se_logit_mmap(
    search_parameters_t* const search_parameters,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification) {
  // Fetch logit model
  const archive_score_logit_model_t* const classify_model =
      archive_score_matches_se_get_model(search_parameters);
  // Compute probability
  const double pr = archive_score_logit_mmap(
      classify_model,matches_predictors,matches_classification);
  const uint64_t mapq_score = archive_score_probability_to_mapq(pr);
  return (uint8_t)MIN(mapq_score,59);
}
/*
 * SE Scoring Models
 */
void archive_score_matches_se_default(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  // Unmapped
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches==0) return;
  // Local matches
  match_trace_t* const primary_map = matches_get_primary_match(matches);
  if (primary_map->type == match_type_local && primary_map->match_alignment.effective_length < 40) {
    primary_map->mapq_score = 1;
    return;
  }
  // Classify
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  matches_predictors_t matches_predictors;
  matches_predictors_compute_se(&matches_predictors,matches);
  // Score
  const int64_t delta = matches->classification.delta_group;
  const int64_t wdelta_group = matches->classification.wdelta_group;
  if (delta == 0) {
    primary_map->mapq_score = 0;
  } else if (delta == -1) {
    if (wdelta_group > 1) {
      primary_map->mapq_score = 60;
    } else {
      primary_map->mapq_score = archive_score_matches_se_logit_unique(
          search_parameters,&matches_predictors,&matches->classification);
    }
  } else {
    if (wdelta_group > 1) {
      primary_map->mapq_score = 60;
    } else {
      primary_map->mapq_score = archive_score_matches_se_logit_mmap(
          search_parameters,&matches_predictors,&matches->classification);
    }
  }
}
void archive_score_matches_se_classify(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  // Unmapped
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches==0) return;
  // Score Local matches
  match_trace_t* const primary_map = matches_get_primary_match(matches);
  if (primary_map->type == match_type_local) {
    primary_map->mapq_score = 1;
    return;
  }
  // Classify
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  const archive_score_logit_model_t* const classify_model = archive_score_matches_se_get_model(search_parameters);
  matches_predictors_t matches_predictors;
  matches_predictors_compute_se(&matches_predictors,matches);
  // Score
  const int64_t delta = matches->classification.delta_group;
  const int64_t wdelta_group = matches->classification.wdelta_group;
  if (delta == 0) {
    // Ties
    matches_metrics_t* const metrics = &matches->metrics;
    if (metrics->min1_event_distance == metrics->min2_event_distance &&
        metrics->min1_edit_distance != metrics->min2_edit_distance &&
        metrics->max1_swg_score != metrics->max2_swg_score) {
      primary_map->mapq_score = 2;
    } else {
      primary_map->mapq_score = 0;
    }
  } else if (delta == -1) {
    if (wdelta_group > 1) {
      if (wdelta_group >= 10) {
        primary_map->mapq_score = 250;
      } else {
        primary_map->mapq_score = 240 + wdelta_group;
      }
    } else {
      const double pr = archive_score_logit_unique(
          classify_model,&matches_predictors,&matches->classification);
      const uint64_t mapq_score = 110 + round(-10.0*log10(1.0-pr));
      primary_map->mapq_score = (mapq_score < 170) ? mapq_score : 170;
    }
  } else {
    if (wdelta_group > 1) {
      if (wdelta_group >= 10) {
        primary_map->mapq_score = 230;
      } else {
        primary_map->mapq_score = 220 + wdelta_group;
      }
    } else {
      const double pr = archive_score_logit_mmap(
          classify_model,&matches_predictors,&matches->classification);
      const uint64_t mapq_score = 10 + round(-10.0*log10(1.0-pr));
      primary_map->mapq_score = (mapq_score < 70) ? mapq_score : 70;
    }
  }
}
/*
 * Archive Scoring SE
 */
void archive_score_matches_se(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  // Check no-matches/no-alignModel
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  if (matches_get_num_match_traces(matches)==0) return;
  if (search_parameters->match_alignment_model==match_alignment_model_none) return;
  /*
   * Select scoring model
   */
  PROFILE_START(GP_ARCHIVE_SCORE_SE_MATCHES,PROFILE_LEVEL);
  switch (search_parameters->mapq_model_se) {
    case mapq_model_none: break;
    case mapq_model_classify:
      archive_score_matches_se_classify(archive_search,matches);
      break;
    case mapq_model_gem:
      archive_score_matches_se_default(archive_search,matches);
      break;
    case mapq_model_dump_predictors: {
      // Score default
      archive_score_matches_se_default(archive_search,matches);
      // Compute predictors
      matches_predictors_t matches_predictors;
      matches_predictors_compute_se(&matches_predictors,matches);
      // Dump predictors
      char* const sequence_tag = sequence_get_tag(archive_search->sequence);
      matches_predictors_se_print(stdout,sequence_tag,&matches_predictors,&matches->classification);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_SE_MATCHES,PROFILE_LEVEL);
}
