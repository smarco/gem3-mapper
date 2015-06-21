/*
 * PROJECT: GEMMapper
 * FILE: archive_score.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_score.h"
#include "matches_classify.h"

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
GEM_INLINE uint8_t archive_score_scale_probability(
    const double probability,const double sum_probability,
    const uint8_t floor,const uint8_t ceil) {
  double mapq_pr = -10. * log10(1.-(probability/sum_probability));
  if (mapq_pr <= 0.0) return floor;
  if (mapq_pr >= 60.0) return ceil;
  return floor + (uint8_t)((mapq_pr*(double)(ceil-floor))/60.);
}
GEM_INLINE int64_t archive_score_matches_gem_se_cases(
    archive_search_t* const archive_search,matches_t* const matches) {
  /*
   * Score Scale
   *
   *   TODO
   */
  // Parameters
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  const uint64_t read_length = sequence_get_length(&archive_search->sequence);
  const swg_penalties_t* const swg_penalties = &archive_search->as_parameters.search_parameters->swg_penalties;
  const uint64_t max_region_length = archive_search_get_max_region_length(archive_search);
  const uint64_t proper_length = fm_index_get_proper_length(archive_search->archive->fm_index);
  // Sort
  matches_sort_by_distance(matches);
  // Compute metrics
  matches_classify_compute_mvalues(matches,swg_penalties,read_length,max_region_length,proper_length,UINT64_MAX);
  // Zero subdominant matches
  double pr;
  uint64_t i;
  for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  if (matches->max_complete_stratum <= 1) return 0; // Remove ambiguous
  // Classify Ties
  matches_tie_type tie_type;
  if ((tie_type=matches_classify_is_tie(matches))) {
    switch (tie_type) {
      case matches_tie_swg_score:
      case matches_tie_edit_distance:
        return 0;
      case matches_tie_event_distance_delta0:
      case matches_tie_event_distance_delta1:
        pr = matches_classify_delta1(matches);
        return (pr >= 0.98) ? archive_score_scale_probability(pr-0.98,0.02,10,29) : 0;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Classify unique matches
  if (matches_classify_is_unique(matches)) {
    pr = matches_classify_unique(matches);
    return (pr >= 0.998) ? archive_score_scale_probability(pr-0.998,0.002,50,60) : 0;
  }
  // Classify multimaps
  pr = matches_classify_mmaps(matches);
  return (pr >= 0.99) ? archive_score_scale_probability(pr-0.99,0.01,30,49) : 0;
}
GEM_INLINE void archive_score_matches_gem_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Parameters
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Classify
  uint64_t i = 0;
  match[0].mapq_score = archive_score_matches_gem_se_cases(archive_search,matches);
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
  const uint64_t read_length = sequence_get_length(&archive_search->sequence);
  const swg_penalties_t* const swg_penalties = &archive_search->as_parameters.search_parameters->swg_penalties;
  const uint64_t max_region_length = archive_search_get_max_region_length(archive_search);
  const uint64_t proper_length = fm_index_get_proper_length(archive_search->archive->fm_index);
  // Compute metrics
  matches_classify_compute_mvalues(matches,swg_penalties,read_length,max_region_length,proper_length,UINT64_MAX);
  // Zero subdominant matches
  double pr;
  uint64_t i;
  for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  // Ties
  matches_tie_type tie_type;
  if ((tie_type=matches_classify_is_tie(matches))) {
    if (matches->max_complete_stratum <= 1) {
      match[0].mapq_score = 1;
    } else {
      switch (tie_type) {
        case matches_tie_swg_score:
          match[0].mapq_score = 2;
          break;
        case matches_tie_edit_distance:
          match[0].mapq_score = 3;
          break;
        case matches_tie_event_distance_delta0:
        case matches_tie_event_distance_delta1:
          pr = matches_classify_delta1(matches);
          match[0].mapq_score = (pr >= 0.90) ? 10 + (uint8_t)((pr-0.90)*500.0) : 9; // 10-60
          break;
        default:
          GEM_INVALID_CASE();
          break;
      }
    }
  } else if (matches_classify_is_unique(matches)) {
    // Classify unique matches
    if (matches->max_complete_stratum <= 1) {
      match[0].mapq_score = 198;
    } else {
      pr = matches_classify_unique(matches);
      match[0].mapq_score = (pr >= 0.95) ? 200 + (uint8_t)((pr-0.95)*1000.0) : 199; // 200-250
    }
  } else {
    // Classify multimaps
    if (matches->max_complete_stratum <= 1) {
      match[0].mapq_score = 138;
    } else {
      pr = matches_classify_mmaps(matches);
      match[0].mapq_score = (pr >= 0.90) ? 140 + (uint8_t)((pr-0.90)*500.0) : 139; // 140-190
    }
  }
  matches_classify_metrics_print(sequence_get_tag(&archive_search->sequence),matches);
}
GEM_INLINE void archive_score_matches_logit_values_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  GEM_NOT_IMPLEMENTED();
}
/*
 * SE Scoring
 */
GEM_INLINE void archive_score_matches_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Check number of matches
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches==0) return;
  PROF_START(GP_ARCHIVE_SCORE_SE_MATCHES);
  // Select scoring model
  switch (archive_search->select_parameters->mapq_model) {
    case mapq_model_none:
      break;
    case mapq_model_gem:
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
  // Check number of matches
  const uint64_t num_matches = vector_get_used(paired_matches->matches);
  if (num_matches==0) return;
  PROF_START(GP_ARCHIVE_SCORE_PE_MATCHES);
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
