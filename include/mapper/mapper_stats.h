/*
 * PROJECT: GEMMapper
 * FILE: mapper_stats.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef MAPPER_STATS_H_
#define MAPPER_STATS_H_

#include "utils/essentials.h"
#include "system/profiler_timer.h"

/*
 * Archive Search Stats
 */
typedef struct {
  /* PE */
  gem_counter_t unique_template_size; // Distribution of the observed template size (for uniq matches)
  /* Performance Stats */
  gem_timer_t generate_candidates_timer;
  uint64_t generate_candidates_total_bases;
  gem_timer_t extend_candidates_timer;
  uint64_t extend_candidates_total_candidates;
  uint64_t extend_candidates_total_bases;
} mapper_stats_t;

/*
 * Setup
 */
mapper_stats_t* mapper_stats_new();
void mapper_stats_clear(mapper_stats_t* const mapper_stats);
void mapper_stats_delete(mapper_stats_t* const mapper_stats);

/*
 * Template Length
 */
void mapper_stats_template_init(
    mapper_stats_t* const search_stats,
    const uint64_t template_length_min,
    const uint64_t template_length_max);
void mapper_stats_template_length_sample(
    mapper_stats_t* const search_stats,
    const uint64_t template_length);

uint64_t mapper_stats_template_length_get_num_samples(mapper_stats_t* const search_stats);
uint64_t mapper_stats_template_length_get_ci_min_samples(
    mapper_stats_t* const search_stats,
    const uint64_t margin_error);
bool mapper_stats_template_length_is_reliable(mapper_stats_t* const search_stats);

double mapper_stats_template_length_get_mean(mapper_stats_t* const search_stats);
double mapper_stats_template_length_get_stddev(mapper_stats_t* const search_stats);
uint64_t mapper_stats_template_length_get_expected_max(mapper_stats_t* const search_stats);
uint64_t mapper_stats_template_length_get_expected_min(mapper_stats_t* const search_stats);
double mapper_stats_template_length_get_sigma_dev(
    mapper_stats_t* const search_stats,
    const uint64_t template_length);


#endif /* MAPPER_STATS_H_ */
