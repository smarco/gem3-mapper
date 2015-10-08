/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_control.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering_control.h"
#include "archive_score.h"
#include "matches_classify.h"

/*
 * Control delta error adjustment
 */
GEM_INLINE void approximate_search_adjust_max_differences_using_strata(
    approximate_search_t* const search,matches_t* const matches) {
  const uint64_t max_differences = search->max_complete_error;
  const uint64_t delta = search->as_parameters->complete_strata_after_best_nominal;
  /*
   * If delta parameter is set (and is below the maximum number of mismatches),
   * finds the minimum non zero stratum (mnzs) and adjusts
   * the maximum number of mismatches to (mnzs+delta)
   */
  if (delta < max_differences) {
    const int64_t fms = matches_metrics_get_min_edit_distance(&matches->metrics);
    if (fms>=0 && fms+delta < max_differences) {
      search->max_complete_error = fms+delta;
    }
  }
}
/*
 * Compute predictors wrapper
 */
GEM_INLINE void asearch_compute_predictors(
    approximate_search_t* const search,matches_t* const matches,
    matches_predictors_t* const predictors) {
  const uint64_t read_length = search->pattern.key_length;
  const swg_penalties_t* const swg_penalties = &search->as_parameters->search_parameters->swg_penalties;
  const uint64_t max_region_length = search->region_profile.max_region_length;
  const uint64_t proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const uint64_t max_complete_stratum = search->max_complete_stratum;
  const uint64_t num_zero_regions = search->region_profile.num_zero_regions;
  matches_classify_compute_predictors(matches,predictors,swg_penalties,
      read_length,max_region_length,proper_length,max_complete_stratum,num_zero_regions);
}
/*
 * Control Fulfilled-search
 */
GEM_INLINE bool asearch_fulfilled(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  if (matches==NULL) return false;
  // Determines when the search is done following the mapping criteria
  switch (search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast: {
      if (matches->max_complete_stratum <= 1) return false;
      const matches_class_t matches_class = matches_classify(matches);
      switch (matches_class) {
        case matches_class_tie_swg_score:
        case matches_class_tie_edit_distance:
        case matches_class_tie_event_distance:
          return (matches->metrics.min1_edit_distance <= 1);
        case matches_class_unique: {
          matches_predictors_t predictors;
          asearch_compute_predictors(search,matches,&predictors);
          return (matches_classify_unique(&predictors) >= MATCHES_UNIQUE_CI);
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
/*
 * Control Filter-ahead
 */
GEM_INLINE bool asearch_filter_ahead_candidates(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const search_parameters = actual_parameters->search_parameters;
  // Determines when the search is done following the mapping criteria
  switch (search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
      return true;
    case mapping_adaptive_filtering_thorough:
      return false;
    case mapping_adaptive_filtering_complete:
      return actual_parameters->complete_strata_after_best_nominal < search->max_complete_error;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
/*
 * Control DFA States
 */
GEM_INLINE bool asearch_control_trigger_boost(approximate_search_t* const search,matches_t* const matches) {
  const matches_class_t matches_class = matches_classify(matches);
  if (matches_class==matches_class_unmapped) {
    return true;
  } else if (matches_class==matches_class_unique) {
    const uint64_t num_zero_regions = search->region_profile.num_zero_regions;
    const bool min_depth = search->max_complete_stratum > num_zero_regions+1;
    const uint64_t proper_length = search->archive->fm_index->proper_length;
    const uint64_t max_region_length = search->region_profile.max_region_length;
    if ((!min_depth && max_region_length > 2*proper_length) || (max_region_length > 3*proper_length)) {
      return true;
    }
    return false;
  } else {
    return false;
  }
}
/*
 * Control Filtering-Thorough
 *   Best efficiency, keeping good-quality results (search for a "good" match)
 */
GEM_INLINE void asearch_control_fast_next_state(
    approximate_search_t* const search,const approximate_search_state_t processing_step,matches_t* const matches) {
  switch (processing_step) {
    case asearch_exact_filtering_adaptive:
      if (search->search_state==asearch_exact_matches) {
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,1);
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,1);
        return;
      }
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      // Boost?
      if (search->search_state==asearch_no_regions) {
        search->search_state = asearch_exact_filtering_boost;
        return;
      }
      search->search_state = asearch_control_trigger_boost(search,matches) ? asearch_exact_filtering_boost : asearch_end;
      break;
    case asearch_exact_filtering_boost: {
      // Unbounded alignment
      search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
      if (matches_is_mapped(matches) || search_parameters->unbounded_alignment==unbounded_alignment_never) {
        search->search_state = asearch_end;
      } else {
        search->search_state = asearch_unbounded_alignment;
      }
      break;
    }
    case asearch_unbounded_alignment:
      search->search_state = asearch_end;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Control Filtering-Thorough
 *   Explore thoroughly all regions/candidates. Spends more computing resources
 *   as to guarantee better-quality results (increasing likelihood of mappings)
 */
GEM_INLINE void asearch_control_thorough_next_state(
    approximate_search_t* const search,
    const approximate_search_state_t processing_step,matches_t* const matches) {
  switch (processing_step) {
    case asearch_exact_filtering_adaptive:
      if (search->search_state==asearch_exact_matches) {
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,1);
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,1);
        return;
      }
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,(matches!=NULL && matches_is_mapped(matches))?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      // Boost
      if (search->search_state==asearch_no_regions) {
        search->search_state = asearch_exact_filtering_boost;
        return;
      }
      search->search_state = asearch_control_trigger_boost(search,matches) ? asearch_exact_filtering_boost : asearch_end;
      break;
    case asearch_exact_filtering_boost: {
      // Unbounded alignment
      search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
      if (matches_is_mapped(matches) || search_parameters->unbounded_alignment==unbounded_alignment_never) {
        search->search_state = asearch_end;
      } else {
        search->search_state = asearch_unbounded_alignment;
      }
      break;
    }
    case asearch_unbounded_alignment:
      search->search_state = asearch_end;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Control Complete
 */
GEM_INLINE void asearch_control_complete_next_state(
    approximate_search_t* const search,const approximate_search_state_t processing_step,
    matches_t* const matches) {
  GEM_NOT_IMPLEMENTED();
}
/*
 * Control HUB
 */
GEM_INLINE void asearch_control_next_state(
    approximate_search_t* const search,const approximate_search_state_t processing_step,
    matches_t* const matches) {
  // Determines when the search is done following the mapping criteria
  switch (search->as_parameters->search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
      asearch_control_fast_next_state(search,processing_step,matches);
      break;
    case mapping_adaptive_filtering_thorough:
      asearch_control_thorough_next_state(search,processing_step,matches);
      break;
    case mapping_adaptive_filtering_complete:
      asearch_control_complete_next_state(search,processing_step,matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
