/*
 * PROJECT: GEMMapper
 * FILE: matches_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches/matches_classify.h"
#include "matches/matches_cigar.h"

/*
 * Classify
 */
void matches_classify(matches_t* const matches) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  // match_trace_t* const match = matches_get_match_trace_buffer(matches);
  // Classify
  if (metrics->total_matches_sampled == 0) {
    matches->matches_class = matches_class_unmapped;
  } else if (metrics->total_matches_sampled == 1) {
    matches->matches_class = matches_class_unique;
  } else if (metrics->max1_swg_score == metrics->max2_swg_score ||
             metrics->min1_edit_distance == metrics->min2_edit_distance ||
             metrics->min1_counter_value == metrics->min2_counter_value) {
          // matches_cigar_cmp(matches->cigar_vector,match,matches->cigar_vector,match+1)==0
    matches->matches_class = matches_class_tie_d0;
  } else if (metrics->min1_edit_distance+1 == metrics->min2_edit_distance ||
             metrics->min1_counter_value+1 == metrics->min2_counter_value) {
    matches->matches_class = matches_class_tie_d1;
  } else {
    matches->matches_class = matches_class_mmap;
  }
}
void paired_matches_classify(paired_matches_t* const paired_matches) {
  // Parameters
  matches_metrics_t* const metrics = &paired_matches->metrics;
  // Unmapped
  const uint64_t num_paired_maps = paired_matches_get_num_maps(paired_matches);
  if (num_paired_maps == 0) {
    paired_matches->paired_matches_class = paired_matches_class_unmapped;
    return;
  }
  // Check both ends MAPQ-quality
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
  match_trace_t* const primary_end1 = paired_map_get_match_end1(paired_matches,paired_map);
  match_trace_t* const primary_end2 = paired_map_get_match_end2(paired_matches,paired_map);
  // Unique
  if (num_paired_maps == 1) {
    paired_matches->paired_matches_class = paired_matches_class_unique;
    return;
  }
  // Tie-d0
  if (metrics->max1_swg_score == metrics->max2_swg_score ||
      metrics->min1_edit_distance == metrics->min2_edit_distance ||
      metrics->min1_counter_value == metrics->min2_counter_value) {
    paired_matches->paired_matches_class = paired_matches_class_tie_d0;
    return;
  }
  match_trace_t* const subdominant_end1 = paired_map_get_match_end1(paired_matches,paired_map+1);
  match_trace_t* const subdominant_end2 = paired_map_get_match_end2(paired_matches,paired_map+1);
  if (matches_cigar_cmp(matches_end1->cigar_vector,primary_end1,matches_end1->cigar_vector,subdominant_end1)==0 &&
      matches_cigar_cmp(matches_end2->cigar_vector,primary_end2,matches_end2->cigar_vector,subdominant_end2)==0) {
    paired_matches->paired_matches_class = paired_matches_class_tie_d0;
    GEM_INVALID_CASE();
    return;
  }
  // Tie-d1
  if (metrics->min1_edit_distance+1 == metrics->min2_edit_distance ||
      metrics->min1_counter_value+1 == metrics->min2_counter_value) {
    paired_matches->paired_matches_class = paired_matches_class_tie_d1;
    return;
  }
  // MultiMap
  paired_matches->paired_matches_class = paired_matches_class_mmap;
}
