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
uint8_t archive_score_matches_se_default_ties(
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  // Classify ties
  const double pr = matches_classify_logit_ties(
      &logit_model_single_end_default,matches_predictors,match_predictors);
  if (pr < MATCHES_MIN_CI) return 0;
  if (pr < MATCHES_TIES_CI) return 1;
  return archive_score_probability_scale(pr-MATCHES_TIES_CI,1.-MATCHES_TIES_CI,3,29);
}
uint8_t archive_score_matches_se_default_mmap(
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  // Classify multimaps
  const double pr = matches_classify_logit_mmaps(
      &logit_model_single_end_default,matches_predictors,match_predictors);
  if (pr < MATCHES_MIN_CI) return 0;
  if (pr < MATCHES_MMAPS_CI) return 2;
  return archive_score_probability_scale(pr-MATCHES_MMAPS_CI,1.-MATCHES_MMAPS_CI,30,49);
}
uint8_t archive_score_matches_se_default_unique(
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  const double pr = matches_classify_logit_unique(
      &logit_model_single_end_default,matches_predictors,match_predictors);
  if (pr < MATCHES_MIN_CI) return 0;
  if (pr < MATCHES_UNIQUE_CI) return 2;
  return archive_score_probability_scale(pr-MATCHES_UNIQUE_CI,1.-MATCHES_UNIQUE_CI,50,60);
}
/*
 * SE Scoring Models
 */
void archive_score_matches_se_stratify(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  // Unmapped
  if (matches->matches_class==matches_class_unmapped) {
    return;
  }
  // Score subdominant matches (MAPQ=0)
  match_trace_t* const match = matches_get_match_trace_buffer(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  uint64_t i;
  for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  // Local
  if (match[0].type == match_type_local) {
    match[0].mapq_score = 1;
    return;
  }
  // Classify
  matches_predictors_t matches_predictors;
  match_predictors_t match_predictors;
  matches_predictors_compute_se(&matches_predictors,matches);
  match_predictors_compute_se(&match_predictors,matches,match);
  switch (matches->matches_class) {
    case matches_class_unmapped:
      break;
    case matches_class_tie_d0:
      match[0].mapq_score = 2;
      break;
    case matches_class_tie_d1: {
      const double pr = matches_classify_logit_ties(
          &logit_model_single_end_default,&matches_predictors,&match_predictors);
      match[0].mapq_score = (pr >= 0.90) ? 80 + (uint8_t)((pr-0.90)*500.0) : 79;
      break;
    }
    case matches_class_mmap: {
      const double pr = matches_classify_logit_mmaps(
          &logit_model_single_end_default,&matches_predictors,&match_predictors);
      match[0].mapq_score = (pr >= 0.90) ? 140 + (uint8_t)((pr-0.90)*500.0) : 139;
      break;
    }
    case matches_class_unique: {
      const double pr = matches_classify_logit_unique(
          &logit_model_single_end_default,&matches_predictors,&match_predictors);
      match[0].mapq_score = (pr >= 0.90) ? 200 + (uint8_t)((pr-0.95)*1000.0) : 199;
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void archive_score_matches_se_default(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  /*
   * Classify
   *   50-60   Unique
   *   30-49   MMaps
   *    2-29   Ties (Distinguishable)
   *       1   Low quality
   *       0   Uncertain & subdominant
   */
  // Unmapped
  if (matches->matches_class == matches_class_unmapped) {
    return;
  }
  // Score subdominant matches (MAPQ=0)
  match_trace_t* const match = matches_get_match_trace_buffer(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  uint64_t i;
  for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  // Classify
  matches_predictors_t matches_predictors;
  match_predictors_t match_predictors;
  switch (matches->matches_class) {
    case matches_class_tie_d0:
      break;
    case matches_class_tie_d1:
      matches_predictors_compute_se(&matches_predictors,matches);
      match_predictors_compute_se(&match_predictors,matches,match);
      match[0].mapq_score = archive_score_matches_se_default_ties(&matches_predictors,&match_predictors);
      break;
    case matches_class_mmap:
      matches_predictors_compute_se(&matches_predictors,matches);
      match_predictors_compute_se(&match_predictors,matches,match);
      match[0].mapq_score = archive_score_matches_se_default_mmap(&matches_predictors,&match_predictors);
      break;
    case matches_class_unique:
      matches_predictors_compute_se(&matches_predictors,matches);
      match_predictors_compute_se(&match_predictors,matches,match);
      match[0].mapq_score = archive_score_matches_se_default_unique(&matches_predictors,&match_predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  matches_metrics_set_mapq(&matches->metrics,match[0].mapq_score);
//  if ( //matches->matches_class==matches_class_unmapped ||
//     ( // matches->matches_class!=matches_class_tie_d0 &&
//      matches->matches_class==matches_class_tie_d1 &&
//      match[0].mapq_score == 0)) {
//    sequence_print(stdout,&archive_search->sequence);
//  }
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
  // Sort by distance (whichever it's selected)
  matches_sort_by_distance(matches);
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
      match_predictors_t match_predictors;
      match_trace_t* const match = matches_is_mapped(matches) ?
          matches_get_match_trace_buffer(matches) : NULL;
      matches_predictors_compute_se(&matches_predictors,matches);
      match_predictors_compute_se(&match_predictors,matches,match);
      // Dump predictors
      matches_predictors_se_print(
          stdout,sequence_get_tag(&archive_search->sequence),
          matches->matches_class,&matches_predictors,&match_predictors);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_SE_MATCHES,PROFILE_LEVEL);
}
