/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_control.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_control.h"
#include "matches/matches_classify.h"
#include "matches/matches_classify_logit.h"
#include "matches/matches_classify_logit_models.h"

/*
 * Search Limits
 */
void asearch_control_adjust_max_differences_using_strata(
    approximate_search_t* const search,
    matches_t* const matches) {
  const uint64_t max_differences = search->max_complete_error;
  const uint64_t delta = search->search_parameters->complete_strata_after_best_nominal;
  /*
   * Control delta error adjustment
   *   If delta parameter is set (and is below the maximum number of mismatches),
   *   finds the minimum non zero stratum (mnzs) and adjusts
   *   the maximum number of mismatches to (mnzs+delta)
   */
  if (delta < max_differences) {
    const int64_t fms = matches_metrics_get_min_edit_distance(&matches->metrics);
    if (fms>=0 && fms+delta < max_differences) {
      search->max_complete_error = fms+delta;
    }
  }
}
/*
 * Search Control Stages
 */
bool asearch_control_filter_ahead_candidates(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  const search_parameters_t* const search_parameters = search->search_parameters;
  // Determines when the search is done following the mapping criteria
  switch (search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
      return true;
    case mapping_adaptive_filtering_thorough:
      return false;
    case mapping_adaptive_filtering_complete:
      return search_parameters->complete_strata_after_best_nominal < search->max_complete_error;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
/*
 * Search Fulfilled & Predictors
 */
void asearch_control_compute_predictors(
    approximate_search_t* const search,
    matches_t* const matches,
    matches_predictors_t* const predictors) {
  matches_predictors_compute(matches,predictors,&search->metrics,search->max_complete_stratum);
}
bool asearch_control_fulfilled(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  const search_parameters_t* const search_parameters = search->search_parameters;
  if (matches==NULL) return false;
  // Determines when the search is done following the mapping criteria
  switch (search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast: {
      if (matches->max_complete_stratum <= 1) return false;
      const matches_class_t matches_class = matches_classify(matches);
      switch (matches_class) {
        case matches_class_tie_d0:
        case matches_class_tie_d1:
          return (matches->metrics.min1_edit_distance <= 1);
        case matches_class_unique: {
          matches_predictors_t predictors;
          asearch_control_compute_predictors(search,matches,&predictors);
          const double probability =
              matches_classify_logit_unique(&predictors,&logit_model_single_end_default);
          return (probability >= MATCHES_UNIQUE_CI);
        }
        case matches_class_unmapped:
        case matches_class_mmap:
        default:
          return false;
      }
      return false;
    }
    case mapping_adaptive_filtering_complete:
      return search->max_complete_stratum > search->max_complete_error;
    case mapping_adaptive_filtering_thorough:
    case mapping_fixed_filtering_complete:
      return false;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}

