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
uint8_t archive_score_probability_to_mapq(
    const double probability,
    const double sum_probability) {
  const double mapq = -10. * log10(1.-(probability/sum_probability));
  if (mapq > 60.) {
    return 60;
  } else if (mapq < 0.) {
    return 0;
  } else {
    return (uint8_t) mapq;
  }
}
uint8_t archive_score_probability_scale(
    const double probability,
    const double sum_probability,
    const uint8_t floor,
    const uint8_t ceil) {
//  double mapq_pr = -10. * log10(1.-(probability/sum_probability));
//  if (mapq_pr <= 0.0) return floor;
//  if (mapq_pr >= 60.0) return ceil;
//  return floor + (uint8_t)((mapq_pr*(double)(ceil-floor))/60.);
  const double range = ceil-floor;
  return floor + (uint8_t)round((range*probability)/sum_probability);
}
/*
 * SE Score Categories
 */
uint8_t archive_score_matches_se_default_ties(matches_predictors_t* const predictors) {
  // Classify ties
  const double pr = matches_classify_logit_ties(predictors,&logit_model_single_end_default);
  if (pr < MATCHES_MIN_CI) return 0;
  if (pr < MATCHES_TIES_CI) return 1;
  return archive_score_probability_scale(pr-MATCHES_TIES_CI,1.-MATCHES_TIES_CI,2,29);
}
uint8_t archive_score_matches_se_default_mmap(matches_predictors_t* const predictors) {
  // Classify multimaps
  const double pr = matches_classify_logit_mmaps(predictors,&logit_model_single_end_default);
  if (pr < MATCHES_MIN_CI) return 0;
  if (pr < MATCHES_MMAPS_CI) return 1;
  return archive_score_probability_scale(pr-MATCHES_MMAPS_CI,1.-MATCHES_MMAPS_CI,30,49);
}
uint8_t archive_score_matches_se_default_unique(matches_predictors_t* const predictors) {
  const double pr = matches_classify_logit_unique(predictors,&logit_model_single_end_default);
  if (pr < MATCHES_MIN_CI) return 0;
  if (pr < MATCHES_UNIQUE_CI) return 1;
  return archive_score_probability_scale(pr-MATCHES_UNIQUE_CI,1.-MATCHES_UNIQUE_CI,50,60);
}
/*
 * SE Scoring Models
 */
void archive_score_matches_se_stratify(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  /*
   * Classify
   *   200-250   Unique (high likelihood)
   *   199       Unique (low likelihood)
   *
   *   140-190   MMaps (high likelihood)
   *   139       MMaps (low likelihood)
   *
   *   80-130    Tie Event (high likelihood)
   *   79        Tie Event (low likelihood)
   *
   *   4         Tie Edit
   *   3         Tie SWG
   *   2         Tie Indistinguishable
   *   1         Ambiguous (not enough resolution mcs <= 1)
   *   0         Subdominant
   */
  // Parameters
  matches_predictors_t predictors;
  match_trace_t* const match = matches_get_match_trace_buffer(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  uint64_t i;
  // Score subdominant matches (MAPQ=0)
  for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  // Compute TP prob
  if (match[0].type == match_type_local) {
    match[0].mapq_score = 1;
  } else {
    switch (matches->matches_class) {
      case matches_class_unmapped:
        break;
      case matches_class_tie_d0:
        match[0].mapq_score = 2;
        break;
      case matches_class_tie_d1: {
        archive_search_se_compute_predictors(archive_search,matches,&predictors);
        const double pr = matches_classify_logit_ties(&predictors,&logit_model_single_end_default);
        match[0].mapq_score = (pr >= 0.90) ? 80 + (uint8_t)((pr-0.90)*500.0) : 79;
        break;
      }
      case matches_class_mmap: {
        archive_search_se_compute_predictors(archive_search,matches,&predictors);
        const double pr = matches_classify_logit_mmaps(&predictors,&logit_model_single_end_default);
        match[0].mapq_score = (pr >= 0.90) ? 140 + (uint8_t)((pr-0.90)*500.0) : 139;
        break;
      }
      case matches_class_unique: {
        archive_search_se_compute_predictors(archive_search,matches,&predictors);
        const double pr = matches_classify_logit_unique(&predictors,&logit_model_single_end_default);
        match[0].mapq_score = (pr >= 0.90) ? 200 + (uint8_t)((pr-0.95)*1000.0) : 199;
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
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
  // Parameters
  matches_predictors_t predictors;
  match_trace_t* const match = matches_get_match_trace_buffer(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  uint64_t i;
  // Score subdominant matches (MAPQ=0)
  for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  // Compute TP prob
  switch (matches->matches_class) {
    case matches_class_unmapped:
    case matches_class_tie_d0:
      break;
    case matches_class_tie_d1:
      archive_search_se_compute_predictors(archive_search,matches,&predictors);
      match[0].mapq_score = archive_score_matches_se_default_ties(&predictors);
      break;
    case matches_class_mmap:
      archive_search_se_compute_predictors(archive_search,matches,&predictors);
      match[0].mapq_score = archive_score_matches_se_default_mmap(&predictors);
      break;
    case matches_class_unique:
      archive_search_se_compute_predictors(archive_search,matches,&predictors);
      match[0].mapq_score = archive_score_matches_se_default_unique(&predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  matches_metrics_set_mapq(&matches->metrics,match[0].mapq_score);
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
  matches->matches_class = matches_classify(matches);
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
      matches_predictors_t predictors;
      archive_search_se_compute_predictors(archive_search,matches,&predictors);
      // Dump predictors
      matches_predictors_se_print(
          stdout,sequence_get_tag(&archive_search->sequence),
          matches->matches_class,&predictors);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_SE_MATCHES,PROFILE_LEVEL);
}
