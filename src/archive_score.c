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
GEM_INLINE void archive_score_matches_se_sm(archive_search_t* const archive_search,matches_t* const matches) {
  // Switch alignment model
  const search_parameters_t* const search_parameters = archive_search->search_actual_parameters.search_parameters;
  if (search_parameters->alignment_model == alignment_model_gap_affine) {
    // Assume sorted descending SWG-score // TODO Sort by score
    const uint64_t num_matches = vector_get_used(matches->global_matches);
    match_trace_t* const match = vector_get_mem(matches->global_matches,match_trace_t);
    // Allocate MAPQ
    // TODO Weight with sequence-quality scores
    mm_stack_t* const mm_stack = archive_search->mm_stack;
    mm_stack_push_state(mm_stack);
    double* const r = mm_stack_calloc(mm_stack,num_matches,double,false);
    // TODO Impl. expected best map to lower score of noisy maps
    const double reference_score = match[0].swg_score;
    double accum_mapq = reference_score;
    const double a = 10.0;
    const double b = (double)-search_parameters->swg_penalties.generic_mismatch_score;
    r[0] = reference_score;
    uint64_t i;
    for (i=1;i<num_matches;++i) {
      double div = pow(a,(double)(reference_score-match[i].swg_score)/b);
      r[i] = (double)match[i-1].swg_score/div;
      accum_mapq += r[i];
    }
    // Calculate final MAPQ scores
    for (i=0;i<num_matches;++i) {
      const double mapq = -10. * log10(1.-(r[i]/accum_mapq));
      if (mapq > 60.) {
        match[i].mapq_score = 60;
      } else if (mapq < 0.) {
        match[i].mapq_score = 0;
      } else {
        match[i].mapq_score = mapq;
      }
    }
    // Free
    mm_stack_pop_state(mm_stack,false);
  } else {
    GEM_NOT_IMPLEMENTED();
  }
}
GEM_INLINE void archive_score_matches_pe_sm(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Switch alignment model
  const search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  if (search_parameters->alignment_model == alignment_model_gap_affine) {
    matches_t* const matches_end1 = paired_matches_end1(paired_matches);
    matches_t* const matches_end2 = paired_matches_end2(paired_matches);
    // Score end1
    archive_score_matches_se_sm(archive_search_end1,matches_end1);
    // Score end2
    archive_score_matches_se_sm(archive_search_end1,matches_end2);
    // Score pairs
    const uint64_t num_matches = vector_get_used(paired_matches->concordant_matches);
    paired_match_t* const paired_match = vector_get_mem(paired_matches->concordant_matches,paired_match_t);
    uint64_t i;
    // COUNTER_PRINT(stderr,&paired_matches->unique_template_size,NULL,"bp",true);
    for (i=0;i<num_matches;++i) { // TODOq
//      const uint64_t num_samples = COUNTER_GET_NUM_SAMPLES(&paired_matches->unique_template_size);
//      const uint64_t mean = (num_samples > 1000) ?
//          : (search_parameters->max_template_length-search_parameters->min_template_length)/2;
      paired_match[i].mapq_score =
          (paired_match[i].match_end1->mapq_score+paired_match[i].match_end2->mapq_score)/2;
    }
  } else {
    GEM_NOT_IMPLEMENTED();
  }
}
/*
 * PR Scoring
 */
GEM_INLINE void archive_score_matches_se_pr(archive_search_t* const archive_search,matches_t* const matches) {
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
GEM_INLINE void archive_score_matches_pe_pr(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  GEM_NOT_IMPLEMENTED();
}
/*
 * Li Scoring
 */
GEM_INLINE void archive_score_matches_se_li(archive_search_t* const archive_search,matches_t* const matches) {
  GEM_NOT_IMPLEMENTED();
}
GEM_INLINE void archive_score_matches_pe_li(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  GEM_NOT_IMPLEMENTED();
}
/*
 * Heath Scoring
 */
GEM_INLINE void archive_score_matches_se_heath(archive_search_t* const archive_search,matches_t* const matches) {
  GEM_NOT_IMPLEMENTED();
}
GEM_INLINE void archive_score_matches_pe_heath(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  GEM_NOT_IMPLEMENTED();
}
/*
 * SE Scoring
 */
GEM_INLINE void archive_score_matches_se(
    archive_search_t* const archive_search,matches_t* const matches) {
  switch (archive_search->select_parameters->mapq_model) {
    case mapq_model_none:
      break;
    case mapq_model_li:
      archive_score_matches_se_li(archive_search,matches);
      break;
    case mapq_model_heath:
      archive_score_matches_se_heath(archive_search,matches);
      break;
    case mapq_model_sm:
      archive_score_matches_se_sm(archive_search,matches);
      break;
    case mapq_model_pr:
      archive_score_matches_se_pr(archive_search,matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * PE Scoring
 */
GEM_INLINE void archive_score_matches_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  switch (archive_search_end1->select_parameters->mapq_model) {
    case mapq_model_none:
      break;
    case mapq_model_li:
      archive_score_matches_pe_li(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_heath:
      archive_score_matches_pe_heath(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_sm:
      archive_score_matches_pe_sm(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_pr:
      archive_score_matches_pe_pr(archive_search_end1,archive_search_end2,paired_matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
