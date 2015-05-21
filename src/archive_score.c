/*
 * PROJECT: GEMMapper
 * FILE: archive_score.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_score.h"

/*
 * Score Cases
 */
GEM_INLINE int64_t archive_score_case(const uint64_t best_distance,const uint64_t mcs,const uint64_t total_matches) {
  /*
   * Noise
   */
  if (best_distance > mcs && best_distance >= 10) return 12;
  if (total_matches > 1) {
    if (total_matches==2) {
      if (mcs >= 4) return 26;
      if (mcs == 2) return 9;
      if (mcs <= 1) return 11;
      if (best_distance >= 4) return 10;
    } else {
      return 0;
    }
  }
  if (total_matches == 1) {
    if (best_distance>=10) return 8;
    if (best_distance>=6 && mcs<=2) return 7;
    if (best_distance==2 && mcs==2) return 6;
    if (best_distance==1 && mcs==2) return 5;
    if (best_distance==1 && mcs==1) return 4;
    if (best_distance==0 && mcs==1) return 3;
    if (mcs == 1) return 2;
  }
  /*
   * Signal
   */
  if (total_matches == 1) {
    if (best_distance < mcs) {
      const uint64_t diff = mcs - (best_distance+1);
      if (diff > 0) return 60;
      if (diff == 0 && mcs >= 3) return 59;
    }
    if (best_distance >= mcs) {
      if (mcs >= 4 && best_distance < 10) return 58;
    }
  }
  // Return unknown
  return -1;
}
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
//  return 100.0 * r_value/accumulated_r;
}
/*
 * GEM (Exponential Relative Score)
 */
GEM_INLINE void archive_score_matches_gem_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Matches Parameters
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  /*
   * Case Scoring (Remove Noise)
   */
  matches_sort_by_distance(matches);
  const uint64_t num_counters = vector_get_used(matches->counters);
  uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  const uint64_t mcs = matches->max_complete_stratum;
  uint64_t i = 0, total_sub = 0;
  while (i < num_counters) total_sub += counters[i++];
  /*
   * Case Scoring (Add signal)
   */
  if (total_sub == 1) {
    const uint64_t best_distance = match[0].distance;
    if (best_distance < mcs) {
      const uint64_t diff = mcs - (best_distance+1);
      if (diff > 0) {
        match[0].mapq_score = 60;
        for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
        return;
      }
      if (diff == 0 && mcs >= 3) {
        match[0].mapq_score = 59;
        for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
        return;
      }
    }
    if (best_distance >= mcs) {
      if (mcs >= 4 && best_distance < 10) {
        match[0].mapq_score = 58;
        for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
        return;
      }
    }
  }
  /*
   * Difference Exponential Score
   */
  const uint64_t read_length = sequence_get_length(&archive_search->sequence);
  const search_parameters_t* const search_parameters = archive_search->search_actual_parameters.search_parameters;
  const int32_t gap_open_score = search_parameters->swg_penalties.gap_open_score;
  const int32_t gap_extension_score = search_parameters->swg_penalties.gap_extension_score;
  const int32_t generic_mismatch_score = search_parameters->swg_penalties.generic_mismatch_score;
  const int32_t generic_match_score = search_parameters->swg_penalties.generic_match_score;
  const int32_t perfect_score = read_length*generic_match_score;
  const double exponential_factor = 1.5;
  // Compute R values
  mm_stack_t* const mm_stack = archive_search->mm_stack;
  mm_stack_push_state(mm_stack);
  double* const r = mm_stack_calloc(mm_stack,num_matches,double,false);
  double accum_r = 0.0;
  // Compute the sum of the R values
  int32_t max_distance = 0, min_distance = INT32_MAX;
  for (i=0;i<num_matches;++i) {
    r[i] = archive_score_diff_exponential(perfect_score,match[i].swg_score,exponential_factor);
    accum_r += r[i];
    if (match[i].distance > max_distance) max_distance = match[i].distance;
    if (match[i].distance < min_distance) min_distance = match[i].distance;
  }
  // Account non-decoded matches (1+1230:9923 || 0:0:1:0+900:10000)
  uint64_t num_counter_matches = 0;
  if (max_distance+1<num_counters) {
    for (i=max_distance+1;i<num_counters;++i) {
      const int32_t average_penality = (gap_open_score + i*gap_extension_score + i*generic_mismatch_score) / 2;
      const int32_t subdominant_score = perfect_score + average_penality;
      const double subdominant_r = archive_score_diff_exponential(perfect_score,subdominant_score,exponential_factor);
      accum_r += counters[i] * subdominant_r;
      num_counter_matches += counters[i];
    }
  }
  // Manual noise-curation (1+0 || 0:0:0:1+0) || (0:0+1:0 || 0:0+0:0:0:0:1)
  if (num_matches > 0) {
    const uint64_t mcs = matches->max_complete_stratum;
    if (num_matches == 1 && mcs == 1) { // (1+0)
      const int32_t undiscovered_distance = max_distance + 1;
      const int32_t undiscovered_avg_penality =
          (gap_open_score + undiscovered_distance*gap_extension_score + undiscovered_distance*generic_mismatch_score) / 2;
      const int32_t undiscovered_score = perfect_score + undiscovered_avg_penality;
      const double undiscovered_r = archive_score_diff_exponential(perfect_score,undiscovered_score,exponential_factor);
      accum_r += undiscovered_r; // Account for the ghost match
    }
    if (min_distance >= mcs || mcs-min_distance==1) { // (0:0+1:0 || 0:0+0:0:0:0:1)
      const int32_t undiscovered_distance = mcs;
      const int32_t undiscovered_avg_penality =
          (gap_open_score + undiscovered_distance*gap_extension_score + undiscovered_distance*generic_mismatch_score) / 2;
      const int32_t undiscovered_score = perfect_score + undiscovered_avg_penality;
      const double undiscovered_r = archive_score_diff_exponential(perfect_score,undiscovered_score,exponential_factor);
      accum_r += undiscovered_r; // Account for the ghost match
    }
  }
  // Calculate final MAPQ scores
  match[0].mapq_score = archive_score_probability_to_mapq(r[0],accum_r);
  for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
  // Free
  mm_stack_pop_state(mm_stack,false);
}
GEM_INLINE void archive_score_matches_gem_case_se(archive_search_t* const archive_search,matches_t* const matches) {
  // Matches Parameters
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  /*
   * Case Scoring (Remove Noise)
   */
  matches_sort_by_distance(matches);
  const uint64_t num_counters = vector_get_used(matches->counters);
  uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  const uint64_t mcs = matches->max_complete_stratum;
  uint64_t i = 0, total_matches = 0;
  while (i < num_counters) total_matches += counters[i++];
  const int64_t score_case = archive_score_case(match[0].distance,mcs,total_matches);
  if (score_case != -1) {
    match[0].mapq_score = score_case;
    for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
    return;
  } else {
    for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
  }
}
GEM_INLINE void archive_score_matches_gem_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  const search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  const uint64_t num_paired_matches = vector_get_used(paired_matches->matches);
  paired_match_t* const paired_match = vector_get_mem(paired_matches->matches,paired_match_t);
  uint64_t i;
  if (num_paired_matches == 0) {
    archive_score_matches_gem_se(archive_search_end1,paired_matches->matches_end1);
    archive_score_matches_gem_se(archive_search_end2,paired_matches->matches_end2);
    return;
  }
  // Multimaps
  if (num_paired_matches > 1) {
    for (i=0;i<num_paired_matches;++i) paired_match[i].mapq_score = 0.0;
    return;
  }
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
///*
// * GEM Stratifying Cases
// */
//GEM_INLINE void archive_score_matches_gem_case_se(archive_search_t* const archive_search,matches_t* const matches) {
//  // Parameters
//  match_trace_t* const match = matches_get_match_traces(matches);
//  const uint64_t num_matches = matches_get_num_match_traces(matches);
//  // Calculate final MAPQ scores
//  if (num_matches > 0) {
//    // Nullify all
//    uint64_t i;
//    for (i=0;i<num_matches;++i) match[i].mapq_score = 0;
//    // MAPQ=0 (multi-maps with the same SWG-score)
//    if (num_matches > 1 && match[0].swg_score==match[1].swg_score) {
//      match[0].mapq_score = 0;
//      return;
//    }
//    // MAPQ=1-9 (multi-maps in the same stratum)
//    uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
//    const uint64_t num_counters = vector_get_used(matches->counters);
//    uint64_t distance = 0;
//    while (distance < num_counters) {
//      if (counters[distance]!=0) {
//        if (counters[distance] > 1) {
//          const uint64_t count = counters[match[0].distance]-2;
//          match[0].mapq_score = count>=8 ? 9 : 1+count;
//          return;
//        }
//        ++distance;
//        break;
//      }
//      ++distance;
//    }
//
//    // MAPQ=10-19 (multi-maps in the (+i) strata)
//    matches_sort_by_distance(matches);
//
//    // Stratify
//    uint64_t total_sub = 0;
//    while (distance < num_counters && counters[distance]==0) ++distance;
//    const uint64_t sub_dis = distance;
//    while (distance < num_counters) total_sub += counters[distance++];
//    if (total_sub > 0) {
////      const uint64_t mcs = matches->max_complete_stratum;
//      const uint64_t diff = sub_dis - match[0].distance;
//      const uint64_t diff_10 = (diff>=10) ? 9 : diff;
//      uint64_t mapq_score;
//      mapq_score = 60 + 10*total_sub + diff_10;
//      match[0].mapq_score = (mapq_score >= 255) ? 254 : mapq_score;
//      return;
//    }
////    // Barrier
////    total_sub = 0; distance = 0;
////    while (distance < num_counters) {
////      total_sub += counters[distance];
////      if (total_sub > 1) {
////        match[0].mapq_score = (distance >= 9) ? 19 : 10 + distance;
////        return;
////      }
////      ++distance;
////    }
//
//    const uint64_t mcs = matches->max_complete_stratum;
//    // MAPQ=20-39 (1 match, mapped beyond the complete strata)
//    if (match[0].distance >= mcs) {
//      if (mcs >= 4 && match[0].distance < 10) {
//        match[0].mapq_score = 59; return;
//      }
//      if (mcs >= 3 && 3 <= match[0].distance && match[0].distance <= 5) {
//        match[0].mapq_score = 57; return;
//      }
//      match[0].mapq_score = (mcs >= 19) ? 20 : 39 - mcs; // Great!!
//      // Stratify
////      const uint64_t f = (mcs>=5) ? 4 : mcs;
////      const uint64_t mapq = 60 + 5*match[0].distance + f;
////      match[0].mapq_score = (mapq > 255) ? 254 : mapq;
//      return;
//    }
//    // MAPQ=40-59 (1 match, mapped into the complete strata)
//    if (match[0].distance < mcs) {
//      const uint64_t diff = mcs - (match[0].distance+1);
//      if (diff > 0) {
//        match[0].mapq_score = 60; return;
//      }
//      if (diff == 0 && mcs >= 3) {
//        match[0].mapq_score = 58; return;
//      }
//      match[0].mapq_score = (diff >= 16) ? 56 : 40 + diff;
////      // Stratify
////      const uint64_t f = (diff>=5) ? 4 : diff;
////      match[0].mapq_score = 60 + 5*mcs + f;
//      return;
//    }
//    // MAPQ=60 (the rest)
//    match[0].mapq_score = 60;
//  }
//}
///*
// * BWA (like)
// */
//GEM_INLINE void archive_score_matches_bwa_se(archive_search_t* const archive_search,matches_t* const matches) {
//  // Parameters
//  match_trace_t* const match = matches_get_match_traces(matches);
//  const uint64_t num_matches = matches_get_num_match_traces(matches);
//  // Calculate final MAPQ scores
//  if (num_matches > 0) {
//    // Read
//    const uint64_t read_length = sequence_get_length(&archive_search->sequence);
//    // Scores
//    const search_parameters_t* const search_parameters = archive_search->search_actual_parameters.search_parameters;
//    // const int32_t generic_mismatch_score = -search_parameters->swg_penalties.generic_mismatch_score;
//    const double generic_match_score = search_parameters->swg_penalties.generic_match_score;
//    // Score the primary one
//    const int32_t sub_swg = (num_matches > 1) ? match[1].swg_score : 0;
//    const double identity = (double)(read_length - match[0].distance)/(double)read_length;
//    double adj = (read_length < 0.5) ? 1. : 3 / log(read_length);
//    adj *= identity * identity;
//    int32_t mapq = (int32_t) (6.02 * (double)(match[0].swg_score - sub_swg) / generic_match_score * adj * adj + 0.499);
//    // printf("mapq = %d adj = %f %f\n",mapq,adj,3 / log(read_length));
//    // Adjust MAPQ using sub-optimal matches
//    uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
//    const uint64_t num_counters = vector_get_used(matches->counters);
//    const uint64_t max_counters = MIN(num_counters,match[0].distance+3);
//    uint64_t i, sub_count = 0;
//    for (i=0;i<max_counters;++i) sub_count += counters[i];
//    if (sub_count > 1) mapq -= (int32_t)(4.343 * log(sub_count+1) + .499);
//    // Set MAPQ
//    if (mapq > 60) mapq = 60;
//    if (mapq < 0) mapq = 0;
//    match[0].mapq_score = mapq;
//    // Nullify the rest
//    for (i=1;i<num_matches;++i) {
//      match[i].mapq_score = 0;
//    }
//  }
//}
/*
 * Logistic Regression
 */
GEM_INLINE void archive_score_matches_logit(archive_search_t* const archive_search,matches_t* const matches) {
  // Matches Parameters
  match_trace_t* const match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  matches_sort_by_distance(matches);
  const uint64_t num_counters = vector_get_used(matches->counters);
  uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  const uint64_t mcs = matches->max_complete_stratum;
  uint64_t i = 0, fs_matches = 0, sub_matches = 0, sub_distance = 0;
  while (i < num_counters && counters[i]==0) i++;
  if (i < num_counters) {
    fs_matches = counters[i++];
    // Get to next non-zero stratum
    while (i < num_counters && counters[i]==0) i++;
    if (i < num_counters) sub_distance = i-match[0].distance;
    while (i < num_counters) sub_matches += counters[i++];
//    // Distance
//    fprintf(stdout,"%lu\t",match[0].distance);
//    // SWG
//    fprintf(stdout,"%d\t",match[0].swg_score);
//    // MCS
//    fprintf(stdout,"%lu\t",mcs);
//    // matches_FS
//    fprintf(stdout,"%lu\t",fs_matches);
//    // sub_matches
//    fprintf(stdout,"%lu\t",sub_matches);
//    // sub_distance
//    fprintf(stdout,"%lu\n",sub_distance);
    const double lr_factor =  -4.56139 +
        match[0].distance * 0.08973 +
        match[0].swg_score * 0.11009 +
        mcs * 0.46871 +
        fs_matches * -2.92813 +
        sub_matches * -0.08716 +
        sub_distance * -0.03892;
    const double pr = 1.0 / (1.0 + (1/exp(lr_factor)));
    match[0].mapq_score = pr*100.0; //archive_score_probability_to_mapq(pr,1.0);
    // Nullify the rest
    for (i=1;i<num_matches;++i) match[i].mapq_score = 0;
  }
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
  const search_parameters_t* const search_parameters = archive_search->search_actual_parameters.search_parameters;
  if (search_parameters->alignment_model != alignment_model_gap_affine) { return; /* GEM_NOT_IMPLEMENTED(); */ }
  // Select scoring model
  switch (archive_search->select_parameters->mapq_model) {
    case mapq_model_none: break;
    case mapq_model_gem:
      archive_score_matches_gem_se(archive_search,matches);
      break;
    case mapq_model_gem_case:
      archive_score_matches_gem_case_se(archive_search,matches);
      break;
    case mapq_model_logit:
      archive_score_matches_logit(archive_search,matches);
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
  const search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  if (search_parameters->alignment_model != alignment_model_gap_affine) { return; /* GEM_NOT_IMPLEMENTED(); */ }
  // Select scoring model
  switch (archive_search_end1->select_parameters->mapq_model) {
    case mapq_model_none: break;
    case mapq_model_gem:
      archive_score_matches_gem_pe(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_gem_case:
    case mapq_model_logit:
      GEM_NOT_IMPLEMENTED();
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SCORE_PE_MATCHES);
}
