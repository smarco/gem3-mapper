/*
 * PROJECT: GEMMapper
 * FILE: archive_score.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_score.h"

/*
 * SM Scoring Utils
 */
GEM_INLINE double archive_score_diff_exponential(
    const int32_t reference_score,const int32_t match_score,double exp_coefficient) {
  return (double)reference_score*exp((double)(match_score-reference_score)*exp_coefficient);
}
GEM_INLINE uint8_t archive_score_probability_to_mapq(const double probability,const double sum_probability) {
   const double mapq = -10. * log10(1.-(probability/sum_probability));
  if (mapq > 60.) {
    return 60;
  } else if (mapq < 0.) {
    return 0;
  } else {
    return (uint8_t) mapq;
  }
}
///*
// * Score Cases
// */
//GEM_INLINE int64_t archive_score_mapq_unique(matches_t* const matches) {
//  const double pr = matches_classify_unique(matches,matches->max_complete_stratum);
//  if (pr >= 0.98) return 60;
//  // Return unknown
//  return -1;
//}
//GEM_INLINE int64_t archive_score_mapq_ambiguous(matches_t* const matches) {
//  // Isolate Pure-Noise
//  const double pr = matches_classify_ambiguous(matches,matches->max_complete_stratum);
//  if (pr <= 0.20) return 0;
//  // Return unknown
//  return -1;
//}
/*
 * GEM Score
 */
//GEM_INLINE double archive_score_matches_compute_r_noise(
//    matches_t* const matches,swg_penalties_t* const swg_penalties,
//    const double exponential_factor,const uint64_t read_length) {
//  // SWG weights
//  const int32_t gap_open_score = swg_penalties->gap_open_score;
//  const int32_t gap_extension_score = swg_penalties->gap_extension_score;
//  const int32_t generic_mismatch_score = swg_penalties->generic_mismatch_score;
//  const int32_t generic_match_score = swg_penalties->generic_match_score;
//  const int32_t perfect_score = read_length*generic_match_score;
//  // Counters
//  const int32_t max_distance = matches->max_counter_value;
//  const int32_t min_distance = matches->min_counter_value;
//  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
//  const uint64_t num_counters = vector_get_used(matches->counters);
//  // Account non-decoded matches (1+1230:9923 || 0:0:1:0+900:10000)
//  double accum_r = 0.0;
//  uint64_t i;
//  if (max_distance+1 < num_counters) {
//    for (i=max_distance+1;i<num_counters;++i) {
//      const int32_t average_penality = (gap_open_score + i*gap_extension_score + i*generic_mismatch_score) / 2;
//      const int32_t subdominant_score = perfect_score + average_penality;
//      const double subdominant_r = archive_score_diff_exponential(perfect_score,subdominant_score,exponential_factor);
//      accum_r += counters[i] * subdominant_r;
//    }
//  }
//  // Manual noise-curation (1+0 || 0:0:0:1+0) || (0:0+1:0 || 0:0+0:0:0:0:1)
//  const uint64_t num_matches = matches_counters_get_total_count(matches);
//  if (num_matches > 0) {
//    const uint64_t mcs = matches->max_complete_stratum;
//    if (num_matches == 1 && mcs == 1) { // (1+0)
//      const int32_t undiscovered_distance = max_distance + 1;
//      const int32_t undiscovered_avg_penality =
//          (gap_open_score + undiscovered_distance*gap_extension_score + undiscovered_distance*generic_mismatch_score) / 2;
//      const int32_t undiscovered_score = perfect_score + undiscovered_avg_penality;
//      const double undiscovered_r = archive_score_diff_exponential(perfect_score,undiscovered_score,exponential_factor);
//      accum_r += undiscovered_r; // Account for the ghost match
//    }
//    if (min_distance >= mcs || mcs-min_distance==1) { // (0:0+1:0 || 0:0+0:0:0:0:1)
//      const int32_t undiscovered_distance = mcs;
//      const int32_t undiscovered_avg_penality = (gap_open_score + undiscovered_distance*gap_extension_score
//          + undiscovered_distance*generic_mismatch_score) / 2;
//      const int32_t undiscovered_score = perfect_score + undiscovered_avg_penality;
//      const double undiscovered_r = archive_score_diff_exponential(perfect_score,undiscovered_score,exponential_factor);
//      accum_r += undiscovered_r; // Account for the ghost match
//    }
//  }
//  return accum_r;
//}
//GEM_INLINE void archive_score_matches_gem_se(archive_search_t* const archive_search,matches_t* const matches) {
//  // Matches Parameters
//  match_trace_t* const match = matches_get_match_traces(matches);
//  const uint64_t num_matches = matches_get_num_match_traces(matches);
//  // Sort
//  uint64_t i;
//  matches_sort_by_swg_score(matches);
//  int64_t mapq;
//  /*
//   * Case Scoring (Remove noise)
//   */
//  if ((mapq=archive_score_mapq_ambiguous(matches)) != -1) {
//    match[0].mapq_score = mapq;
//    for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
//    return;
//  }
//  /*
//   * Case Scoring (Add signal)
//   */
//  if ((mapq=archive_score_mapq_unique(matches)) != -1) {
//    match[0].mapq_score = mapq;
//    for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
//    return;
//  }
//  /*
//   * Difference Exponential Score
//   */
//  swg_penalties_t* const swg_penalties = &archive_search->as_parameters.search_parameters->swg_penalties;
//  const uint64_t read_length = sequence_get_length(&archive_search->sequence);
//  const int32_t perfect_score = read_length*swg_penalties->generic_match_score;
//  const double exponential_factor = 1.5;
//  // Compute R values
//  const double primary_r = archive_score_diff_exponential(perfect_score,match[0].swg_score,exponential_factor);
//  double accum_r = primary_r;
//  for (i=1;i<num_matches;++i) {
//    accum_r += archive_score_diff_exponential(perfect_score,match[i].swg_score,exponential_factor);
//  }
//  accum_r += archive_score_matches_compute_r_noise(matches,swg_penalties,exponential_factor,read_length);
//  // Calculate final MAPQ scores
//  match[0].mapq_score = archive_score_probability_to_mapq(primary_r,accum_r);
//  if (match[0].mapq_score > 57) match[0].mapq_score = 57;
//  for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
//}
GEM_INLINE int64_t archive_score_matches_gem_se_cases(matches_t* const matches) {
  /*
   * Score Scale
   *
   *   60 - Unique mappings
   *   [10,60) - MMaps range
   *   3 - Ambiguous Mapping (Coudn't be classified with confidence)
   *   2 - Ambiguous Mapping (High confidence that the TP cannot be resolved)
   *   1 - Ties (SWG-ties or EditDistance-ties)
   *   0 - Sub-dominant matches (wrt SWG score)
   */
  // Remove ties
  double pr = matches_classify_ties(matches);
  if (pr == 1.0) return 1;
  // Classify Ambiguous (Remove hard/fuzzy to classify)
  pr = matches_classify_ambiguous(matches,matches->max_complete_stratum);
  if (pr <= 0.98) return 2;
  // Classify unique matches (num_matches==1)
  pr = matches_classify_unique(matches,matches->max_complete_stratum);
  if (pr >= 0.999) return 60;
  // Classify remaining multimaps
  pr = matches_classify_mmaps(matches,matches->max_complete_stratum);
  if (pr >= 0.98) return (int64_t)(round((pr-0.98)*2500.0))+10;
  return 3;
}
GEM_INLINE void archive_score_matches_gem_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Parameters
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Sort
  matches_sort_by_swg_score(matches);
  // matches_metrics_print(matches);
  uint64_t i = 0;
  match[0].mapq_score = archive_score_matches_gem_se_cases(matches);
  for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
}
GEM_INLINE void archive_score_matches_gem_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  const search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  const uint64_t num_paired_matches = vector_get_used(paired_matches->matches);
  paired_match_t* const paired_match = vector_get_mem(paired_matches->matches,paired_match_t);
  uint64_t i;
  if (num_paired_matches == 0) {
    archive_score_matches_gem_se(archive_search_end1,paired_matches->matches_end1);
    archive_score_matches_gem_se(archive_search_end2,paired_matches->matches_end2);
    return;
  }
  // Multimaps
  if (num_paired_matches == 1) {
    paired_match[0].mapq_score = 60;
    for (i=1;i<num_paired_matches;++i) paired_match[i].mapq_score = 0.0;
  } else {
    // Score Parameters
    const uint64_t read_length_1 = sequence_get_length(&archive_search_end1->sequence);
    const uint64_t read_length_2 = sequence_get_length(&archive_search_end2->sequence);
    const uint64_t total_read_length = read_length_1 + read_length_2;
    const int32_t generic_match_score = search_parameters->swg_penalties.generic_match_score;
    const int32_t perfect_score = total_read_length*generic_match_score;
    // Compute R values
    mm_stack_t* const mm_stack = archive_search_end1->mm_stack;
    mm_stack_push_state(mm_stack);
    double* const r = mm_stack_calloc(mm_stack,num_paired_matches,double,false);
    double accum_r = 0.0;
    for (i=0;i<num_paired_matches;++i) {
      if (paired_match[i].pair_orientation==pair_orientation_discordant) {
        r[i] = 0.0;
      } else {
        // Compute differential r-value
        const int32_t match_score = paired_match[i].match_end1->swg_score + paired_match[i].match_end2->swg_score;
        r[i] = archive_score_diff_exponential(perfect_score,match_score,1.0);
        accum_r += r[i];
      }
    }
    // Calculate final MAPQ scores
    if (accum_r==0.0) accum_r = 1.0;
    for (i=0;i<num_paired_matches;++i) {
      paired_match[i].mapq_score = archive_score_probability_to_mapq(r[i],accum_r);
    }
    // Free
    mm_stack_pop_state(mm_stack,false);
  }
}
/*
 * Logit Values
 */
GEM_INLINE void archive_score_matches_logit_values_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Parameters
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Sort
  matches_sort_by_swg_score(matches);
  /*
   * Score Scale
   *
   *   190-250 - Unique mappings
   *   120-180 - MMaps range
   *    50-110 - Ambiguous Mappings
   *         2 - Coudn't be classified (with confidence)
   *         1 - Ties (SWG-ties or EditDistance-ties)
   *         0 - Sub-dominant matches (wrt SWG score)
   */
  // Sub-dominant matches
  uint64_t i;
  for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  // Remove ties
  double pr = matches_classify_ties(matches);
  if (pr == 1.0) {
    match[0].mapq_score = 1;
    matches_metrics_print(matches);
    return;
  }
  // Classify Ambiguous (Remove hard/fuzzy to classify)
  pr = matches_classify_ambiguous(matches,matches->max_complete_stratum);
  if (pr <= 0.98) {
    match[0].mapq_score = 50 + archive_score_probability_to_mapq(pr-0.98,0.02);
    matches_metrics_print(matches);
    return;
  }
  // Classify unique matches (num_matches==1)
  pr = matches_classify_unique(matches,matches->max_complete_stratum);
  if (pr >= 0.999) {
    match[0].mapq_score = 190 + archive_score_probability_to_mapq(pr-0.999,0.001);
    matches_metrics_print(matches);
    return;
  }
  // Classify remaining multimaps
  pr = matches_classify_mmaps(matches,matches->max_complete_stratum);
  if (pr >= 0.98) {
    match[0].mapq_score = 120 + archive_score_probability_to_mapq(pr-0.98,0.001);
    matches_metrics_print(matches);
    return;
  }
  // Otherwise Ambiguous
  match[0].mapq_score = 2;
  matches_metrics_print(matches);
}
GEM_INLINE void archive_score_matches_logit_values_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  GEM_NOT_IMPLEMENTED();
}
/*
 * SE Scoring
 */
GEM_INLINE void archive_score_matches_se(
    archive_search_t* const archive_search,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SCORE_SE_MATCHES);
  // Check number of matches
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches==0) return;
  // Check alignment model
  const search_parameters_t* const search_parameters = archive_search->as_parameters.search_parameters;
  if (search_parameters->alignment_model != alignment_model_gap_affine) { return; /* GEM_NOT_IMPLEMENTED(); */ }
  // Select scoring model
  switch (archive_search->select_parameters->mapq_model) {
    case mapq_model_none: break;
    case mapq_model_gem:
      // archive_score_matches_gem_se(archive_search,matches);
      archive_score_matches_gem_se(archive_search,matches);
      break;
    case mapq_model_logit:
      archive_score_matches_logit_values_se(archive_search,matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SCORE_SE_MATCHES);
}
/*
 * PE Scoring
 */
GEM_INLINE void archive_score_matches_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROF_START(GP_ARCHIVE_SCORE_PE_MATCHES);
  // Check number of matches
  const uint64_t num_matches = vector_get_used(paired_matches->matches);
  if (num_matches==0) return;
  // Check alignment model
  const search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  if (search_parameters->alignment_model != alignment_model_gap_affine) { return; /* GEM_NOT_IMPLEMENTED(); */ }
  // Select scoring model
  switch (archive_search_end1->select_parameters->mapq_model) {
    case mapq_model_none: break;
    case mapq_model_gem:
      archive_score_matches_gem_pe(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_logit:
      archive_score_matches_logit_values_pe(archive_search_end1,archive_search_end2,paired_matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SCORE_PE_MATCHES);
}
