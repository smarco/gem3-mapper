/*
 * PROJECT: GEMMapper
 * FILE: archive_score.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_score.h"

/*
 * SM Scoring
 */
GEM_INLINE double archive_score_diff_exponential(
    const int32_t reference_score,const int32_t match_score,double exp_coefficient) {
  const double e = 2.71828182845904523;
  return (double)reference_score*pow(e,(double)(match_score-reference_score)*exp_coefficient);
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
//  return 100.0 * r_value/accumulated_r;
}
/*
 * SE
 *   // TODO Weight with sequence-quality scores
 */
GEM_INLINE void archive_score_matches_exp_relative_distance_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Parameters
  match_trace_t* const match = vector_get_mem(matches->global_matches,match_trace_t);
  const uint64_t num_matches = vector_get_used(matches->global_matches);
  // Score Parameters
  const uint64_t read_length = sequence_get_length(&archive_search->sequence);
  const uint64_t mcs = matches->max_complete_stratum;
  // Estimate an undiscovered map
  const int32_t perfect_score = read_length;
  const int32_t undiscovered_score = perfect_score - mcs;
  // Compute R values
  mm_stack_t* const mm_stack = archive_search->mm_stack;
  mm_stack_push_state(mm_stack);
  double* const r = mm_stack_calloc(mm_stack,num_matches,double,false);
  const double exponential_factor = 8.0;
  double accum_r = archive_score_diff_exponential(perfect_score,undiscovered_score,exponential_factor);
  int32_t i;
  for (i=0;i<num_matches;++i) {
    r[i] = archive_score_diff_exponential(perfect_score,perfect_score-match[i].distance,exponential_factor);
    accum_r += r[i];
  }
  // Calculate final MAPQ scores
  for (i=0;i<num_matches;++i) {
    match[i].mapq_score = archive_score_probability_to_mapq(r[i],accum_r);
  }
  // Free
  mm_stack_pop_state(mm_stack,false);
}
GEM_INLINE void archive_score_matches_exp_relative_score_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Parameters
  match_trace_t* const match = vector_get_mem(matches->global_matches,match_trace_t);
  const uint64_t num_matches = vector_get_used(matches->global_matches);
  const search_parameters_t* const search_parameters = archive_search->search_actual_parameters.search_parameters;
  // Score Parameters
  const uint64_t read_length = sequence_get_length(&archive_search->sequence);
//  const uint64_t mcs = matches->max_complete_stratum;
  const int32_t gap_open_score = search_parameters->swg_penalties.gap_open_score;
  const int32_t gap_extension_score = search_parameters->swg_penalties.gap_extension_score;
  const int32_t generic_mismatch_score = search_parameters->swg_penalties.generic_mismatch_score;
  const int32_t generic_match_score = search_parameters->swg_penalties.generic_match_score;
  const int32_t perfect_score = read_length*generic_match_score;
  // Estimate an undiscovered map
  const double exponential_factor = 1.0;
//  const int32_t undiscovered_score = match[0].swg_score + average_penality*((int32_t)mcs-(int32_t)match[0].distance);
//  const double undiscovered_r = archive_score_diff_exponential(perfect_score,undiscovered_score,exponential_factor);
  // Compute R values
  mm_stack_t* const mm_stack = archive_search->mm_stack;
  mm_stack_push_state(mm_stack);
  double* const r = mm_stack_calloc(mm_stack,num_matches,double,false);
  double accum_r = 0.0;
//  // Manual noise-curation (1+0)
//  if (match[0].distance==0 && mcs==1) accum_r = undiscovered_r;
//  // Manual noise-curation (0:0+1:0 || 0:0+0:1 || 0:0+0:0:1)
//  if ((int32_t)match[0].distance-(int32_t)mcs >= 0) accum_r = undiscovered_r;
  int32_t i, max_distance = 0;
  for (i=0;i<num_matches;++i) {
    r[i] = archive_score_diff_exponential(perfect_score,match[i].swg_score,exponential_factor);
    accum_r += r[i];
    if (match[i].distance > max_distance) max_distance = match[i].distance;
  }
  // Manual noise-curation (1+1230:9923 || 0:0:1:0+900:10000)
  const uint64_t num_counters = vector_get_used(matches->counters);
  if (max_distance+1<num_counters) {
    uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
    for (i=max_distance+1;i<num_counters;++i) {
      const int32_t average_penality = (gap_open_score + i*gap_extension_score + i*generic_mismatch_score) / 2;
      const int32_t subdominant_score = perfect_score + average_penality;
      const double subdominant_r = archive_score_diff_exponential(perfect_score,subdominant_score,exponential_factor);
      accum_r += counters[i] * subdominant_r;
    }
  }
  // Calculate final MAPQ scores
  for (i=0;i<num_matches;++i) {
    match[i].mapq_score = archive_score_probability_to_mapq(r[i],accum_r);
  }
  // Free
  mm_stack_pop_state(mm_stack,false);
}
/*
 * PE
 */
GEM_INLINE void archive_score_matches_exp_relative_distance_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
//  // Switch alignment model
//  const search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
//  if (search_parameters->alignment_model == alignment_model_gap_affine) {
//    // Score pairs
//    const uint64_t num_matches = vector_get_used(paired_matches->concordant_matches);
//    paired_match_t* const paired_match = vector_get_mem(paired_matches->concordant_matches,paired_match_t);
//    // Allocate R values
//    mm_stack_t* const mm_stack = archive_search_end1->mm_stack;
//    mm_stack_push_state(mm_stack);
//    double* const r = mm_stack_calloc(mm_stack,num_matches,double,false);
//    // Compute R values
//    const uint64_t read_length_end1 = sequence_get_length(&archive_search_end1->sequence);
//    const uint64_t read_length_end2 = sequence_get_length(&archive_search_end2->sequence);
//    const int32_t reference_swg_score = (read_length_end1+read_length_end2) * search_parameters->swg_penalties.generic_match_score;
//    double accum_r = 0.0;
//    const uint64_t num_samples = COUNTER_GET_NUM_SAMPLES(&paired_matches->unique_template_size);
//    if (num_samples >= 100) {
//      const double mean = COUNTER_GET_MEAN(&paired_matches->unique_template_size);
//      const double std_dev = COUNTER_GET_STDDEV(&paired_matches->unique_template_size);
//      uint64_t i;
//      for (i=0;i<num_matches;++i) {
//        // Compute z-score & CDF
//        const double z_score = fabs((double)paired_match[i].template_length-mean) / std_dev;
//        const double p_template_length = 2.0*standard_normal_CDF(-z_score);
//        // Compute R-value
//        const int64_t swg_score = paired_match[i].match_end1->swg_score + paired_match[i].match_end2->swg_score;
//        r[i] = p_template_length * archive_score_matches_sm_match(reference_swg_score,swg_score);
//        accum_r += r[i];
//      }
//    } else {
//      uint64_t i;
//      for (i=0;i<num_matches;++i) {
//        const int64_t swg_score = paired_match[i].match_end1->swg_score + paired_match[i].match_end2->swg_score;
//        r[i] = archive_score_matches_sm_match(reference_swg_score,swg_score);
//        accum_r += r[i];
//      }
//    }
//    // Calculate final MAPQ scores
//    uint64_t i;
//    for (i=0;i<num_matches;++i) {
//      paired_match[i].mapq_score = archive_score_matches_sm_mapq(r[i],accum_r);
//    }
//    // Free
//    mm_stack_pop_state(mm_stack,false);
//  } else {
//    GEM_NOT_IMPLEMENTED();
//  }
}
/*
 * PR Scoring
 */
GEM_INLINE void archive_score_matches_test_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Matches
  const uint64_t num_matches = vector_get_used(matches->global_matches);
  if (num_matches==0) return;
  match_trace_t* const match = vector_get_mem(matches->global_matches,match_trace_t);
  /*
   * Tr1
   */
//  const uint8_t mapq = (40./sqrt((double)count));
//  match[0].mapq_score = mapq;
//  for (i=1;i<num_matches;++i,++count) {
//    if (reference_distance == match[i].distance) {
//      match[i].mapq_score = mapq;
//    } else {
//      match[i].mapq_score = 0;
//    }
//  }
  /*
   * Tr2
   */
//  if (count==1) {
//    if (num_matches==1) {
//      match[0].mapq_score = 60;
//    } else if (next-1 == reference_distance) {
//      match[0].mapq_score = 0;
//    } else {
//      match[0].mapq_score = 60;
//    }
//  } else {
//    match[0].mapq_score = 0;
//  }
  /*
   * Tr3
   */
//  if (num_matches==1) {
//    match[0].mapq_score = 10 + (match[0].swg_score*(100-match[0].distance))/200;
//  } else {
//    match[0].mapq_score = 0;
//  }

  /*
   * SW is the Smith-Waterman scheme,
   * L the read length,
   * n_matches the total number of matches,
   * n_first the number of matches in the first stratum
   * n_sub the number of matches in the second stratum
   */
  const double l = archive_search->sequence.read.length;
  const double n_matches = vector_get_used(matches->global_matches);
  double n_first = 1.0;
  double n_sub = 0.0;
  // Count same distance matches
  uint64_t i = 1;
  uint64_t s0_distance = match[0].distance, s1_distance = 0;
  if (i<num_matches) {
    for (;i<num_matches;++i) {
      if (s0_distance != match[i].distance) break;
      n_first = n_first + 1.;
    }
    if (i<num_matches) {
      s1_distance = match[i].distance;
      for (;i<num_matches;++i) {
        if (s1_distance != match[i].distance) break;
        n_sub = n_sub + 1.;
      }
    }
  }
  /*
   * Tr4
   */
  for (i=0;i<num_matches;++i) {
    const double SW = match[i].swg_score;
    //  1. (SW / L * 60) / (n_matches * n_matches)
    //match[i].mapq_score =  ((double)SW/l * 60.)/(n_matches*n_matches);
    //  2. (SW / L * 60) / (n_first * (1 + n_sub))
    //match[i].mapq_score =  ((double)SW/l * 60.)/(n_first*(1. + n_sub));
    //  3. (SW / L * 60) / (n_matches * n_first)
    //match[i].mapq_score =  ((double)SW/l * 60.)/(n_matches*n_first);
    //  4. (SW / L * 60) / (n_matches * n_matches * n_matches)
    //match[i].mapq_score =  ((double)SW/l * 60.)/(n_matches*n_matches*n_matches);
    //  5. (SW / L * 60) / (n_matches * n_first * (1 + n_sub))
    //match[i].mapq_score =  ((double)SW/l * 60.)/(n_matches*n_first*(1.+n_sub));
    //  6. (SW * 60) / (L* n_matches * n_matches)
    //match[i].mapq_score = ((double)SW*60.) / (l*n_matches*n_matches);
    //  7. (SW * 60) / (L * n_first * (1 + n_sub))
    //match[i].mapq_score = ((double)SW*60.) / (l*n_first*(1.+n_sub));
    //  8. (SW * 60) / (L * n_matches * n_first)
    //match[i].mapq_score = ((double)SW*60.) / (l*n_matches*n_first);
    //  9. (SW * 60) / (L * n_matches * n_matches * n_matches)
    //match[i].mapq_score = ((double)SW*60.) / (l*n_matches*n_matches*n_matches);
    // 10. (SW * 60) / (L * n_matches * n_first * (1 + n_sub))
     match[i].mapq_score = ((double)SW*60.) / (l*n_matches*n_first*(1.+n_sub));
  }
}
GEM_INLINE void archive_score_matches_test_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  GEM_NOT_IMPLEMENTED();
}
/*
 * SE Scoring
 */
GEM_INLINE void archive_score_matches_se(
    archive_search_t* const archive_search,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SCORE_MATCHES);
  // Check number of matches
  const uint64_t num_matches = vector_get_used(matches->global_matches);
  if (num_matches==0) return;
  // Check alignment model
  const search_parameters_t* const search_parameters = archive_search->search_actual_parameters.search_parameters;
  if (search_parameters->alignment_model != alignment_model_gap_affine) { return; /* GEM_NOT_IMPLEMENTED(); */ }
  // Select scoring model
  switch (archive_search->select_parameters->mapq_model) {
    case mapq_model_none:
      break;
    case mapq_model_exp_relative_distance:
      archive_score_matches_exp_relative_distance_se(archive_search,matches);
      break;
    case mapq_model_exp_relative_score:
      archive_score_matches_exp_relative_score_se(archive_search,matches);
      break;
    case mapq_model_test:
      archive_score_matches_test_se(archive_search,matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SCORE_MATCHES);
}
/*
 * PE Scoring
 */
GEM_INLINE void archive_score_matches_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROF_START(GP_ARCHIVE_SCORE_MATCHES);
  // Check number of matches
  const uint64_t num_matches = vector_get_used(paired_matches->concordant_matches);
  if (num_matches==0) return;
  // Check alignment model
  const search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  if (search_parameters->alignment_model != alignment_model_gap_affine) { return; /* GEM_NOT_IMPLEMENTED(); */ }
  // Select scoring model
  switch (archive_search_end1->select_parameters->mapq_model) {
    case mapq_model_none:
      break;
    case mapq_model_exp_relative_distance:
      //archive_score_matches_exp_relative_distance_pe(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_exp_relative_score:
      //archive_score_matches_exp_relative_score_pe(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_test:
      archive_score_matches_test_pe(archive_search_end1,archive_search_end2,paired_matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SCORE_MATCHES);
}
