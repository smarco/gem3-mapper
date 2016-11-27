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
#include "matches/classify/matches_classify_logit.h"
#include "matches/classify/matches_classify_logit_models.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Scoring Utils
 */
uint8_t archive_score_probability_scale(
    const double probability,
    const double sum_probability,
    const uint8_t floor,
    const uint8_t ceil) {
  const double range = ceil-floor;
  return floor + (uint8_t)round((range*probability)/sum_probability);
}
/*
 * SE Predictors Model
 */
const matches_classify_logit_model_t* archive_score_matches_se_get_model(
    search_parameters_t* const search_parameters) {
  if (search_parameters->mapping_mode==mapping_adaptive_filtering_fast) {
    return &logit_model_single_end_fast;
  } else {
    return &logit_model_single_end_sensitive;
  }
}
/*
 * SE Score Categories
 */
uint8_t archive_score_matches_se_default_mmaps_d1(
    const matches_classify_logit_model_t* const classify_model,
    const matches_predictors_t* const matches_predictors) {
  // Classify ties
  const double pr = matches_classify_logit_mmaps_d1(classify_model,matches_predictors);
  if (pr < MATCHES_MIN_CI) return 4;
  if (pr < MATCHES_TIES_CI) return 7;
  return archive_score_probability_scale(pr-MATCHES_TIES_CI,1.-MATCHES_TIES_CI,10,29);
}
uint8_t archive_score_matches_se_default_mmap(
    const matches_classify_logit_model_t* const classify_model,
    const matches_predictors_t* const matches_predictors) {
  // Classify multimaps
  const double pr = matches_classify_logit_mmaps(classify_model,matches_predictors);
  if (pr < MATCHES_MIN_CI) return 4;
  if (pr < MATCHES_MMAPS_CI) return 8;
  return archive_score_probability_scale(pr-MATCHES_MMAPS_CI,1.-MATCHES_MMAPS_CI,30,49);
}
uint8_t archive_score_matches_se_default_unique(
    const matches_classify_logit_model_t* const classify_model,
    const matches_predictors_t* const matches_predictors) {
  const double pr = matches_classify_logit_unique(classify_model,matches_predictors);
  if (pr < MATCHES_MIN_CI) return 4;
  if (pr < MATCHES_UNIQUE_CI) return 9;
  return archive_score_probability_scale(pr-MATCHES_UNIQUE_CI,1.-MATCHES_UNIQUE_CI,50,60);
}
/*
 * SE Scoring Models
 */
void archive_score_matches_se_default(
    matches_t* const matches,
    const matches_classify_logit_model_t* const classify_model) {
  // Unmapped
  if (matches->matches_class == matches_class_unmapped) return;
  // Score subdominant matches (MAPQ=0)
  match_trace_t* const match = matches_get_primary_match(matches);
  // Local
  if (match[0].type == match_type_local) {
    match[0].mapq_score = 2;
    return;
  }
  // Classify
  matches_predictors_t matches_predictors;
  switch (matches->matches_class) {
    case matches_class_tie_perfect:
      if (matches_get_num_match_traces(matches)==2) {
        match[0].mapq_score = 1;
      }
      break;
    case matches_class_tie:
      match[0].mapq_score = 3;
      break;
    case matches_class_mmap_d1:
      matches_predictors_compute_se(&matches_predictors,matches);
      match[0].mapq_score = archive_score_matches_se_default_mmaps_d1(classify_model,&matches_predictors);
      break;
    case matches_class_mmap:
      matches_predictors_compute_se(&matches_predictors,matches);
      match[0].mapq_score = archive_score_matches_se_default_mmap(classify_model,&matches_predictors);
      break;
    case matches_class_unique:
      matches_predictors_compute_se(&matches_predictors,matches);
      match[0].mapq_score = archive_score_matches_se_default_unique(classify_model,&matches_predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  matches_metrics_set_mapq(&matches->metrics,match[0].mapq_score);
}
void archive_score_matches_se_stratify(
    matches_t* const matches,
    const matches_classify_logit_model_t* const classify_model) {
  // Unmapped
  if (matches->matches_class==matches_class_unmapped) return;
  // Score subdominant matches (MAPQ=0)
  match_trace_t* const match = matches_get_primary_match(matches);
  // Local
  if (match[0].type == match_type_local) {
    match[0].mapq_score = 1;
    return;
  }
  // Classify
  matches_predictors_t matches_predictors;
  matches_predictors_compute_se(&matches_predictors,matches);
  switch (matches->matches_class) {
    case matches_class_unmapped:
      break;
    case matches_class_tie_perfect:
      match[0].mapq_score = 2;
      break;
    case matches_class_tie: {
      match[0].mapq_score = 3;
      break;
    }
    case matches_class_mmap_d1: {
      const double pr = matches_classify_logit_mmaps_d1(classify_model,&matches_predictors);
      if (pr < 0.90) {
        match[0].mapq_score = 39;
      } else {
        match[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,40,90);
      }
      break;
    }
    case matches_class_mmap: {
      const double pr = matches_classify_logit_mmaps(classify_model,&matches_predictors);
      if (pr < 0.90) {
        match[0].mapq_score = 99;
      } else {
        match[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,100,150);
      }
      break;
    }
    case matches_class_unique: {
      const double pr = matches_classify_logit_unique(classify_model,&matches_predictors);
      if (pr < 0.90) {
        match[0].mapq_score = 199;
      } else {
        match[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,200,250);
      }
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
}
#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */
void archive_score_matches_se_profile(matches_t* const matches) {
  // Profile matches classification
  switch (matches->matches_class) {
    case matches_class_tie_perfect: PROF_INC_COUNTER(GT_MATCHES_SE_TIE_PERFECT); break;
    case matches_class_tie: PROF_INC_COUNTER(GT_MATCHES_SE_TIE); break;
    case matches_class_mmap_d1: PROF_INC_COUNTER(GT_MATCHES_SE_MMAP_D1); break;
    case matches_class_mmap: PROF_INC_COUNTER(GT_MATCHES_SE_MMAP); break;
    case matches_class_unique: PROF_INC_COUNTER(GT_MATCHES_SE_UNIQUE); return;
    case matches_class_unmapped: PROF_INC_COUNTER(GT_MATCHES_SE_UNMAPPED); return;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Parameters
  profile_t* const profile = PROF_GET_PROFILE();
  matches_metrics_t* const metrics = &matches->metrics;
  match_trace_t* const primary_match = matches_get_match_traces(matches)[0];
  match_trace_t* const secondary_match = matches_get_match_traces(matches)[1];
  // Compute strata-delta (edit)
  if (primary_match->edit_distance > metrics->min_edit_distance) {
    PROF_INC_COUNTER(GT_MATCHES_SE_EDIT_UNCONSISTENT);
  } else {
    PROF_INC_COUNTER(GT_MATCHES_SE_EDIT_CONSISTENT);
    stats_vector_add(profile->strata_deltas_edit,secondary_match->edit_distance-primary_match->edit_distance,1);
  }
  // Compute strata-delta (swg)
  if (primary_match->swg_score < metrics->max_swg_score) {
    PROF_INC_COUNTER(GT_MATCHES_SE_SWG_UNCONSISTENT);
  } else {
    PROF_INC_COUNTER(GT_MATCHES_SE_SWG_CONSISTENT);
    stats_vector_add(profile->strata_deltas_swg,primary_match->swg_score-secondary_match->swg_score,1);
  }
}
#else
void archive_score_matches_se_profile(matches_t* const matches) {}
#endif
/*
 * Archive Scoring SE
 */
void archive_score_matches_se(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  // Check no-matches
  if (matches_get_num_match_traces(matches)==0) return;
  // Check no-align model
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  if (search_parameters->match_alignment_model==match_alignment_model_none) return;
  // Fetch predictors-model
  const matches_classify_logit_model_t* const classify_model = archive_score_matches_se_get_model(search_parameters);
  // Classify
  matches_classify(matches);
  archive_score_matches_se_profile(matches);
  /*
   * Select scoring model
   */
  PROFILE_START(GP_ARCHIVE_SCORE_SE_MATCHES,PROFILE_LEVEL);
  switch (search_parameters->mapq_model_se) {
    case mapq_model_none: break;
    case mapq_model_classify:
      archive_score_matches_se_stratify(matches,classify_model);
      break;
    case mapq_model_gem:
      archive_score_matches_se_default(matches,classify_model);
      break;
    case mapq_model_dump_predictors: {
      // Score default
      archive_score_matches_se_default(matches,classify_model);
      // Compute predictors
      matches_predictors_t matches_predictors;
      matches_predictors_compute_se(&matches_predictors,matches);
      // Dump predictors
      matches_predictors_se_print(
          stdout,sequence_get_tag(&archive_search->sequence),
          matches->matches_class,&matches_predictors);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_SE_MATCHES,PROFILE_LEVEL);
}
