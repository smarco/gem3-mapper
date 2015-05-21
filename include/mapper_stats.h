/*
 * PROJECT: GEMMapper
 * FILE: mapper_stats.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef MAPPER_STATS_H_
#define MAPPER_STATS_H_

#include "essentials.h"

/*
 * Constants
 */
#define MS_TEMPLATE_LENGTH_DEFAULT_MOE 20

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
GEM_INLINE mapper_stats_t* mapper_stats_new();
GEM_INLINE void mapper_stats_clear(mapper_stats_t* const mapper_stats);
GEM_INLINE void mapper_stats_delete(mapper_stats_t* const mapper_stats);

/*
 * Template Size
 */
GEM_INLINE void mapper_stats_template_length_sample(
    mapper_stats_t* const search_stats,const uint64_t template_length);

GEM_INLINE uint64_t mapper_stats_template_length_get_num_samples(mapper_stats_t* const search_stats);
GEM_INLINE uint64_t mapper_stats_template_length_get_ci_min_samples(
    mapper_stats_t* const search_stats,const uint64_t margin_error);
GEM_INLINE bool mapper_stats_template_length_estimation_within_ci(
    mapper_stats_t* const search_stats,const uint64_t margin_error);

GEM_INLINE double mapper_stats_template_length_get_mean(mapper_stats_t* const search_stats);
GEM_INLINE double mapper_stats_template_length_get_stddev(mapper_stats_t* const search_stats);
GEM_INLINE uint64_t mapper_stats_template_length_get_expected_max(mapper_stats_t* const search_stats);
GEM_INLINE uint64_t mapper_stats_template_length_get_expected_min(mapper_stats_t* const search_stats);


#endif /* MAPPER_STATS_H_ */
