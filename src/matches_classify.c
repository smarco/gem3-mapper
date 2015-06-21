/*
 * PROJECT: GEMMapper
 * FILE: matches_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches_classify.h"

/*
 * Setup
 */
GEM_INLINE void matches_classify_metrics_init(matches_metrics_t* const matches_metrics) {
  // Matches metrics
  matches_metrics->min1_counter_value = UINT32_MAX;
  matches_metrics->min2_counter_value = UINT32_MAX;
  matches_metrics->max_counter_value = 0;
  matches_metrics->min1_edit_distance = UINT32_MAX;
  matches_metrics->min2_edit_distance = UINT32_MAX;
  matches_metrics->max1_swg_score = INT32_MIN;
  matches_metrics->max2_swg_score = INT32_MIN;
}
/*
 * Accessors
 */
GEM_INLINE uint64_t matches_classify_metrics_get_min_distance(matches_t* const matches) {
  return matches->metrics.min1_counter_value;
}
GEM_INLINE uint64_t matches_classify_metrics_get_max_distance(matches_t* const matches) {
  return matches->metrics.max_counter_value;
}
GEM_INLINE uint64_t matches_classify_metrics_get_min_edit_distance(matches_t* const matches) {
  return matches->metrics.min1_edit_distance;
}
GEM_INLINE int32_t matches_classify_metrics_get_max_swg_score(matches_t* const matches) {
  return matches->metrics.max1_swg_score;
}
/*
 * Metrics
 */
GEM_INLINE void matches_classify_metrics_update(
    matches_metrics_t* const matches_metrics,const uint64_t distance,
    const uint64_t edit_distance,const int32_t swg_score) {
  // Max counter
  matches_metrics->max_counter_value = MAX(matches_metrics->max_counter_value,distance);
  // Min counter
  if (distance <= matches_metrics->min1_counter_value) {
    matches_metrics->min2_counter_value = matches_metrics->min1_counter_value;
    matches_metrics->min1_counter_value = distance;
  } else if (distance < matches_metrics->min2_counter_value) {
    matches_metrics->min2_counter_value = distance;
  }
  // Min Edit-distance
  if (edit_distance <= matches_metrics->min1_edit_distance) {
    matches_metrics->min2_edit_distance = matches_metrics->min1_edit_distance;
    matches_metrics->min1_edit_distance = edit_distance;
  } else if (edit_distance < matches_metrics->max2_swg_score) {
    matches_metrics->min2_edit_distance = edit_distance;
  }
  // Max SGW-score
  if (swg_score >= matches_metrics->max1_swg_score) {
    matches_metrics->max2_swg_score = matches_metrics->max1_swg_score;
    matches_metrics->max1_swg_score = swg_score;
  } else if (swg_score > matches_metrics->max2_swg_score) {
    matches_metrics->max2_swg_score = swg_score;
  }
}
GEM_INLINE void matches_classify_metrics_recompute(matches_t* const matches) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  // ReCompute distance metrics
  const match_trace_t* match = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  matches_classify_metrics_init(metrics);
  uint64_t i;
  for (i=0;i<num_matches;++i,++match) {
    matches_classify_metrics_update(metrics,match->distance,match->edit_distance,match->swg_score);
  }
}
GEM_INLINE void matches_classify_compute_mvalues(
    matches_t* const matches,const swg_penalties_t* const swg_penalties,
    const uint64_t read_length,const uint64_t max_region_length,
    uint64_t const proper_length,uint64_t const overriding_mcs) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  // Compute event distances
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches == 0) {
    metrics->first_map_event_distance = read_length;
    metrics->subdominant_event_distance = read_length;
    metrics->first_map_edit_distance = read_length;
    metrics->subdominant_edit_distance = read_length;
    metrics->first_map_swg_score = 0;
    metrics->subdominant_swg_score = 0;
  } else {
    match_trace_t* const match = matches_get_match_traces(matches);
    // Compute primary/subdominant event distance
    metrics->first_map_event_distance = match->distance;
    if (match->distance!=metrics->min1_counter_value) {
      metrics->subdominant_event_distance = metrics->min1_counter_value;
    } else {
      metrics->subdominant_event_distance = metrics->min2_counter_value;
    }
    if (metrics->subdominant_event_distance==UINT32_MAX) {
      metrics->subdominant_event_distance = read_length;
    }
    // Compute primary/subdominant edit distance
    metrics->first_map_edit_distance = match->edit_distance;
    if (match->edit_distance!=metrics->min1_edit_distance) {
      metrics->subdominant_edit_distance = metrics->min1_edit_distance;
    } else {
      metrics->subdominant_edit_distance = metrics->min2_edit_distance;
    }
    if (metrics->subdominant_edit_distance==UINT32_MAX) {
      metrics->subdominant_edit_distance = read_length;
    }
    // Compute primary/subdominant SWG score
    metrics->first_map_swg_score = match->swg_score;
    if (metrics->first_map_swg_score < 0) metrics->first_map_swg_score = 0;
    if (match->swg_score!=metrics->max1_swg_score) {
      metrics->subdominant_swg_score = metrics->max1_swg_score;
    } else {
      metrics->subdominant_swg_score = metrics->max2_swg_score;
    }
    if (metrics->subdominant_swg_score < 0) metrics->subdominant_swg_score = 0;
  }
  // Normalize event distances
  metrics->first_map_event_distance_norm =
      ((double)(read_length-metrics->first_map_event_distance))/((double)read_length);
  metrics->subdominant_event_distance_norm =
      ((double)(read_length-metrics->subdominant_event_distance))/((double)read_length);
  // Normalize edit distances
  metrics->first_map_edit_distance_norm =
      ((double)(read_length-metrics->first_map_edit_distance))/((double)read_length);
  metrics->subdominant_edit_distance_norm =
      ((double)(read_length-metrics->subdominant_edit_distance))/((double)read_length);
  // Normalize SWG scores
  const double swg_norm_factor = swg_penalties->generic_match_score*read_length;
  metrics->first_map_swg_score_norm = (double)metrics->first_map_swg_score/swg_norm_factor;
  metrics->subdominant_swg_score_norm = (double)metrics->subdominant_swg_score/swg_norm_factor;
  // Max region length
  metrics->max_region_length = max_region_length;
  metrics->max_region_length_norm = (double)max_region_length/(double)proper_length;
  // Maximum complete stratum
  metrics->mcs = (overriding_mcs!=UINT64_MAX) ? overriding_mcs : matches->max_complete_stratum;
  // First/Subdominant strata
  metrics->first_stratum_matches =
      (metrics->min1_counter_value==UINT32_MAX) ? 0 : matches_counters_get_count(matches,metrics->min1_counter_value);
  metrics->subdominant_stratum_matches = matches_counters_get_total_count(matches) - metrics->first_stratum_matches;
}
GEM_INLINE matches_tie_type matches_classify_is_tie(matches_t* const matches) {
  matches_metrics_t* const metrics = &matches->metrics;
  if (metrics->max1_swg_score == metrics->max2_swg_score) {
    return matches_tie_swg_score;
  } else if (metrics->min1_edit_distance == metrics->min2_edit_distance) {
    return matches_tie_edit_distance;
  } else if (metrics->first_stratum_matches > 1) {
    return matches_tie_event_distance_delta0;
  } else if (matches->metrics.subdominant_event_distance-matches->metrics.first_map_event_distance == 1) {
    return matches_tie_event_distance_delta1;
  } else {
    return matches_tie_none;
  }
}
GEM_INLINE bool matches_classify_is_unique(matches_t* const matches) {
  matches_metrics_t* const metrics = &matches->metrics;
  return (metrics->first_stratum_matches == 1 && metrics->subdominant_stratum_matches == 0);
}
GEM_INLINE double matches_classify_unique(matches_t* const matches) {
  // Unique: Probability of the first position-match (primary match) of being a true positive
  matches_metrics_t* const metrics = &matches->metrics;
  const double lr_factor = -204.9124 +
      (double)metrics->first_map_event_distance_norm * 187.6882 +
      (double)metrics->first_map_swg_score_norm * 18.1521 +
      (double)metrics->mcs * 5.1507 +
      (double)metrics->max_region_length_norm * -0.6575;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
GEM_INLINE double matches_classify_mmaps(matches_t* const matches) {
  // Classify MMaps wrt the probability of the first position-match (primary match) of being a true positive
  matches_metrics_t* const metrics = &matches->metrics;
  const double lr_factor = -196.3411 +
      (double)metrics->first_map_event_distance_norm * 202.4395 +
      (double)metrics->subdominant_event_distance_norm * -20.9556 +
      (double)metrics->first_map_swg_score_norm * 21.1421 +
      (double)metrics->subdominant_swg_score_norm * -7.9429 +
      (double)metrics->mcs * 3.5875 +
      (double)metrics->subdominant_stratum_matches * 0.0844;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
GEM_INLINE double matches_classify_delta1(matches_t* const matches) {
  // Classify ties wrt the probability of being a true positive
  matches_metrics_t* const metrics = &matches->metrics;
  const double lr_factor = -161.35009 +
      (double)metrics->first_map_edit_distance_norm * 16.98019 +
      (double)metrics->subdominant_edit_distance_norm * -5.63140 +
      (double)metrics->first_map_event_distance_norm * 261.40959 +
      (double)metrics->subdominant_event_distance_norm * -112.04288 +
      (double)metrics->first_map_swg_score_norm * 12.54211 +
      (double)metrics->subdominant_swg_score_norm * -13.87366 +
      (double)metrics->mcs * 1.52146 +
      (double)metrics->max_region_length_norm * 0.30583 +
      (double)metrics->first_stratum_matches * -0.10691 +
      (double)metrics->subdominant_stratum_matches * 0.01256;
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
/*
 * Display
 */
GEM_INLINE void matches_classify_metrics_print(const char* const read_tag,matches_t* const matches) {
  matches_metrics_t* const metrics = &matches->metrics;
  // edit
  fprintf(stdout,"%f\t",metrics->first_map_edit_distance_norm);
  // sub_edit
  fprintf(stdout,"%f\t",metrics->subdominant_edit_distance_norm);
  // event
  fprintf(stdout,"%f\t",metrics->first_map_event_distance_norm);
  // sub_event
  fprintf(stdout,"%f\t",metrics->subdominant_event_distance_norm);
  // swg
  fprintf(stdout,"%f\t",metrics->first_map_swg_score_norm);
  // sub_swg
  fprintf(stdout,"%f\t",metrics->subdominant_swg_score_norm);
  // mcs
  fprintf(stdout,"%lu\t",metrics->mcs);
  // max_region_length
  fprintf(stdout,"%f\t",metrics->max_region_length_norm);
  // fs_matches
  fprintf(stdout,"%lu\t",metrics->first_stratum_matches);
  // sub_matches
  fprintf(stdout,"%lu\t",metrics->subdominant_stratum_matches);
  // delta
  fprintf(stdout,"%lu\t",metrics->subdominant_event_distance-metrics->first_map_event_distance);
  // tag
  fprintf(stdout,"%s\t",read_tag);
  // mapq_score
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches > 0) {
    match_trace_t* const match = matches_get_match_traces(matches);
    fprintf(stdout,"%d\n",match[0].mapq_score);
  } else {
    fprintf(stdout,"\n");
  }
}
