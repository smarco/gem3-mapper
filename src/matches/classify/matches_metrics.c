/*
 * PROJECT: GEMMapper
 * FILE: matches_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches/classify/matches_metrics.h"

/*
 * Setup
 */
void matches_metrics_init(matches_metrics_t* const metrics) {
  // Matches
  metrics->accepted_candidates = 0;
  metrics->accepted_matches = 0;
  metrics->limited_candidates = false;
  // Matches Distance
  metrics->min_event_distance = UINT32_MAX;
  metrics->min_event_distance_count = 0;
  metrics->min_edit_distance = UINT32_MAX;
  metrics->min_edit_distance_count = 0;
  metrics->max_swg_score = INT32_MIN;
  metrics->max_swg_score_count = 0;
  // Matches template-length
  metrics->min_template_length_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
  metrics->min_template_length_sigma_count = 0;
  // MAPQ Score
  metrics->mapq = 0;
}
/*
 * Accessors
 */
uint64_t matches_metrics_get_min_event_distance(matches_metrics_t* const metrics) {
  return metrics->min_event_distance;
}
uint64_t matches_metrics_get_min_edit_distance(matches_metrics_t* const metrics) {
  return metrics->min_edit_distance;
}
int32_t matches_metrics_get_max_swg_score(matches_metrics_t* const metrics) {
  return metrics->max_swg_score;
}
/* */
void matches_metrics_set_proper_length(
    matches_metrics_t* const metrics,
    const double proper_length) {
  metrics->proper_length = proper_length;
}
void matches_metrics_set_read_length(
    matches_metrics_t* const metrics,
    const uint64_t read_length) {
  metrics->read_length = read_length;
}
void matches_metrics_set_swg_match_score(
    matches_metrics_t* const metrics,
    const int32_t swg_match_score) {
  metrics->swg_match_score = swg_match_score;
}
void matches_metrics_set_max_region_length(
    matches_metrics_t* const metrics,
    const uint64_t max_region_length) {
  metrics->max_region_length = max_region_length;
}
void matches_metrics_set_kmer_frequency(
    matches_metrics_t* const metrics,
    const double kmer_frequency) {
  metrics->kmer_frequency = kmer_frequency;
}
void matches_metrics_add_accepted_candidates(
    matches_metrics_t* const metrics,
    const uint64_t accepted_candidates) {
  metrics->accepted_candidates += accepted_candidates;
}
void matches_metrics_set_accepted_candidates(
    matches_metrics_t* const metrics,
    const uint64_t accepted_candidates) {
  metrics->accepted_candidates = accepted_candidates;
}
void matches_metrics_set_limited_candidates(
    matches_metrics_t* const metrics,
    const bool limited_candidates) {
  metrics->limited_candidates = limited_candidates;
}
void matches_metrics_set_mapq(
    matches_metrics_t* const metrics,
    const uint8_t mapq) {
  metrics->mapq = mapq;
}
/*
 * Update
 */
void matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t event_distance,
    const uint64_t edit_distance,
    const int32_t swg_score) {
  // Samples
  ++(matches_metrics->accepted_matches);
  // Min counter
  if (event_distance < matches_metrics->min_event_distance) {
    matches_metrics->min_event_distance = event_distance;
    matches_metrics->min_event_distance_count = 1;
  } else if (event_distance == matches_metrics->min_event_distance) {
    ++(matches_metrics->min_event_distance_count);
  }
  // Min Edit-distance
  if (edit_distance < matches_metrics->min_edit_distance) {
    matches_metrics->min_edit_distance = edit_distance;
    matches_metrics->min_edit_distance_count = 1;
  } else if (edit_distance == matches_metrics->min_edit_distance) {
    ++(matches_metrics->min_edit_distance_count);
  }
  // Max SWG-score
  if (swg_score > matches_metrics->max_swg_score) {
    matches_metrics->max_swg_score = swg_score;
    matches_metrics->max_swg_score_count = 1;
  } else if (swg_score == matches_metrics->max_swg_score) {
    ++(matches_metrics->max_swg_score_count);
  }
}
void paired_matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t event_distance,
    const uint64_t edit_distance,
    const int32_t swg_score,
    const double template_length_sigma) {
  // Update metrics
  matches_metrics_update(matches_metrics,event_distance,edit_distance,swg_score);
  // Min TemplateSize-Sigma
  if (template_length_sigma < matches_metrics->min_template_length_sigma) {
    matches_metrics->min_template_length_sigma = template_length_sigma;
    matches_metrics->min_template_length_sigma_count = 1;
  } else if (template_length_sigma == matches_metrics->min_template_length_sigma) {
    ++(matches_metrics->min_template_length_sigma_count);
  }
}
/*
 * Display
 */
void matches_metrics_print(
    FILE* const stream,
    matches_metrics_t* const matches_metrics) {
  tab_fprintf(stream,"[GEM]>Metrics\n");
  tab_fprintf(stream,"  => Search.Magnitudes\n");
  tab_fprintf(stream,"    => Read.length     %lu\n",matches_metrics->read_length);
  tab_fprintf(stream,"    => Proper.Length   %2.3f\n",matches_metrics->proper_length);
  tab_fprintf(stream,"    => SWG.Match.Score %lu\n",matches_metrics->swg_match_score);
  tab_fprintf(stream,"  => Mappability\n");
  tab_fprintf(stream,"    => Max.Region.length  %lu\n",matches_metrics->max_region_length);
  tab_fprintf(stream,"    => Kmer.frequency     %2.3f\n",matches_metrics->kmer_frequency);
  tab_fprintf(stream,"  => Matches\n");
  tab_fprintf(stream,"    => Accepted.candidates     %lu\n",matches_metrics->accepted_candidates);
  tab_fprintf(stream,"    => Accepted.matches        %lu\n",matches_metrics->accepted_matches);
  tab_fprintf(stream,"    => Min.counter.value       %lu\n",matches_metrics->min_event_distance);
  tab_fprintf(stream,"    => Min.counter.value.count %lu\n",matches_metrics->min_event_distance_count);
  tab_fprintf(stream,"    => Min.edit.distance       %lu\n",matches_metrics->min_edit_distance);
  tab_fprintf(stream,"    => Min.edit.distance.count %lu\n",matches_metrics->min_edit_distance_count);
  tab_fprintf(stream,"    => Max.swg.score           %ld\n",matches_metrics->max_swg_score);
  tab_fprintf(stream,"    => Max.swg.score.count     %lu\n",matches_metrics->max_swg_score_count);
  tab_fprintf(stream,"  => PE.Specific\n");
  tab_fprintf(stream,"    => Min.template.sigma       %2.3f\n",matches_metrics->min_template_length_sigma);
  tab_fprintf(stream,"    => Min.template.sigma.count %lu\n",matches_metrics->min_template_length_sigma_count);
  tab_fprintf(stream,"  => MAPQ\n");
  tab_fprintf(stream,"    => MAPQ %d\n",matches_metrics->mapq);
}
