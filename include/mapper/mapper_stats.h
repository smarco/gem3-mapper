/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Mapper-stats data structure enables keeping track of relevant
 *   magnitudes along the mapping-search (e.g. template length) as to
 *   aid the mapping process to take decision towards accuracy, precision,
 *   performance or simply to report statistics of the process
 */

#ifndef MAPPER_STATS_H_
#define MAPPER_STATS_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"

/*
 * Archive Search Stats
 */
typedef struct {
  /* PE */
  gem_counter_t unique_template_size;          // Distribution of the observed template-size (for uniq matches)
  gem_counter_t unique_template_size_ref;      // Reference template-size to start collecting samples
  /* Performance Stats */
  gem_timer_t generate_candidates_timer;       // TODO
  uint64_t generate_candidates_total_bases;    // TODO
  gem_timer_t extend_candidates_timer;         // TODO
  uint64_t extend_candidates_total_candidates; // TODO
  uint64_t extend_candidates_total_bases;      // TODO
} mapper_stats_t;

/*
 * Setup
 */
mapper_stats_t* mapper_stats_new(void);
void mapper_stats_clear(mapper_stats_t* const mapper_stats);
void mapper_stats_delete(mapper_stats_t* const mapper_stats);

/*
 * Template Length
 */
void mapper_stats_template_init(
    mapper_stats_t* const search_stats,
    const uint64_t template_num_samples,
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
