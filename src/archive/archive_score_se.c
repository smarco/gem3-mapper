/*
 * PROJECT: GEMMapper
 * FILE: archive_score.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_score_se.h"
#include "archive/archive_search_se.h"
#include "matches/matches_classify.h"
#include "matches/matches_classify_logit.h"
#include "matches/matches_classify_logit_models.h"

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
 * SE Score Categories
 */
uint8_t archive_score_matches_se_default_mmaps_d1(const matches_predictors_t* const matches_predictors) {
  // Classify ties
  const double pr = matches_classify_logit_mmaps_d1(
      &logit_model_single_end_default,matches_predictors);
  if (pr < MATCHES_MIN_CI) return 4;
  if (pr < MATCHES_TIES_CI) return 7;
  return archive_score_probability_scale(pr-MATCHES_TIES_CI,1.-MATCHES_TIES_CI,10,29);
}
uint8_t archive_score_matches_se_default_mmap(const matches_predictors_t* const matches_predictors) {
  // Classify multimaps
  const double pr = matches_classify_logit_mmaps(
      &logit_model_single_end_default,matches_predictors);
  if (pr < MATCHES_MIN_CI) return 4;
  if (pr < MATCHES_MMAPS_CI) return 8;
  return archive_score_probability_scale(pr-MATCHES_MMAPS_CI,1.-MATCHES_MMAPS_CI,30,49);
}
uint8_t archive_score_matches_se_default_unique(const matches_predictors_t* const matches_predictors) {
  const double pr = matches_classify_logit_unique(
      &logit_model_single_end_default,matches_predictors);
  if (pr < MATCHES_MIN_CI) return 4;
  if (pr < MATCHES_UNIQUE_CI) return 9;
  return archive_score_probability_scale(pr-MATCHES_UNIQUE_CI,1.-MATCHES_UNIQUE_CI,50,60);
}
/*
 * SE Scoring Models
 */
void archive_score_matches_se_default(
    archive_search_t* const archive_search,
    matches_t* const matches) {
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
      match[0].mapq_score = archive_score_matches_se_default_mmaps_d1(&matches_predictors);
      break;
    case matches_class_mmap:
      matches_predictors_compute_se(&matches_predictors,matches);
      match[0].mapq_score = archive_score_matches_se_default_mmap(&matches_predictors);
      break;
    case matches_class_unique:
      matches_predictors_compute_se(&matches_predictors,matches);
      match[0].mapq_score = archive_score_matches_se_default_unique(&matches_predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  matches_metrics_set_mapq(&matches->metrics,match[0].mapq_score);
}
void archive_score_matches_se_stratify(
    archive_search_t* const archive_search,
    matches_t* const matches) {
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
      const double pr = matches_classify_logit_mmaps_d1(
          &logit_model_single_end_default,&matches_predictors);
      if (pr < 0.90) {
        match[0].mapq_score = 39;
      } else {
        match[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,40,90);
      }
      break;
    }
    case matches_class_mmap: {
      const double pr = matches_classify_logit_mmaps(
          &logit_model_single_end_default,&matches_predictors);
      if (pr < 0.90) {
        match[0].mapq_score = 99;
      } else {
        match[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,100,150);
      }
      break;
    }
    case matches_class_unique: {
      const double pr = matches_classify_logit_unique(
          &logit_model_single_end_default,&matches_predictors);
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
  if (search_parameters->alignment_model==alignment_model_none) return;
  // Classify
  matches_classify(matches);
  /*
   * Select scoring model
   */
  PROFILE_START(GP_ARCHIVE_SCORE_SE_MATCHES,PROFILE_LEVEL);
  switch (search_parameters->mapq_model_se) {
    case mapq_model_none: break;
    case mapq_model_classify:
      archive_score_matches_se_stratify(archive_search,matches);
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
