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

#include "mapper/mapper_stats.h"

/*
 * Debug
 */
#define MAPPING_STATS_LOG false

/*
 * Constants
 */
#define MAPPER_STATS_TEMPLATE_LENGTH_MIN_SAMPLES   100
#define MAPPER_STATS_TEMPLATE_LENGTH_MARGIN_ERROR   20

/*
 * Setup
 */
mapper_stats_t* mapper_stats_new(void) {
  // Alloc
  mapper_stats_t* const mapper_stats = mm_alloc(mapper_stats_t);
  // Init
  mapper_stats_clear(mapper_stats);
  // Return
  return mapper_stats;
}
void mapper_stats_clear(mapper_stats_t* const mapper_stats) {
  // PE
  COUNTER_RESET(&mapper_stats->unique_template_size);
  COUNTER_RESET(&mapper_stats->unique_template_size_ref);
  // Performance Stats
  TIMER_RESET(&mapper_stats->generate_candidates_timer);
  mapper_stats->generate_candidates_total_bases = 0;
  TIMER_RESET(&mapper_stats->extend_candidates_timer);
  mapper_stats->extend_candidates_total_candidates = 0;
  mapper_stats->extend_candidates_total_bases = 0;
}
void mapper_stats_delete(mapper_stats_t* const mapper_stats) {
  mm_free(mapper_stats);
}
/*
 * Template Length
 */
void mapper_stats_template_init(
    mapper_stats_t* const search_stats,
    const uint64_t template_num_samples,
    const uint64_t template_length_min,
    const uint64_t template_length_max) {
  if (template_length_min!=UINT64_MAX &&
      template_length_max!=UINT64_MAX &&
      template_length_min <= template_length_max) {
    gem_counter_t* const tlength = &search_stats->unique_template_size;
    const double mean = ((double)(template_length_min+template_length_max))/2.0;
    const double dev = mean/6.0;
    tlength->samples = template_num_samples;
    tlength->min = template_length_min;
    tlength->max = template_length_max;
    tlength->m_oldM = mean;
    tlength->m_newM = tlength->m_oldM;
    tlength->m_oldS = dev*dev*template_num_samples;
    tlength->m_newS = tlength->m_oldS;
    tlength->total = (uint64_t)((double)template_num_samples * mean);
  }
}
void mapper_stats_template_length_sample(
    mapper_stats_t* const search_stats,
    const uint64_t template_length) {
  COUNTER_ADD(&search_stats->unique_template_size,template_length);
  gem_cond_log(MAPPING_STATS_LOG,
      "[GEM]> MappingStats.templateLength {mean=%f,std_dev=%f} "
      "(samples=%"PRIu64",samples_ci=%"PRIu64")",
      mapper_stats_template_length_get_mean(search_stats),
      mapper_stats_template_length_get_stddev(search_stats),
      mapper_stats_template_length_get_num_samples(search_stats),
      mapper_stats_template_length_get_ci_min_samples(search_stats,20));
}
uint64_t mapper_stats_template_length_get_num_samples(mapper_stats_t* const search_stats) {
  return COUNTER_GET_NUM_SAMPLES(&search_stats->unique_template_size);
}
uint64_t mapper_stats_template_length_get_ci_min_samples(
    mapper_stats_t* const search_stats,
    const uint64_t margin_error) {
  const double moe = margin_error; // Margin error
  const double z = 1.96; // 95%
  const double std_dev = COUNTER_GET_STDDEV(&search_stats->unique_template_size);
  const double factor = (z*std_dev)/moe;
  const uint64_t min_samples = (uint64_t)(factor*factor)+1;
  return min_samples;
}
bool mapper_stats_template_length_is_reliable(mapper_stats_t* const search_stats) {
  const uint64_t num_samples = mapper_stats_template_length_get_num_samples(search_stats);
  const uint64_t num_ci_min_samples =
      mapper_stats_template_length_get_ci_min_samples(search_stats,MAPPER_STATS_TEMPLATE_LENGTH_MARGIN_ERROR);
  return num_samples > num_ci_min_samples && num_samples >= MAPPER_STATS_TEMPLATE_LENGTH_MIN_SAMPLES;
}
double mapper_stats_template_length_get_mean(mapper_stats_t* const search_stats) {
  return COUNTER_GET_MEAN(&search_stats->unique_template_size);
}
double mapper_stats_template_length_get_stddev(mapper_stats_t* const search_stats) {
  return COUNTER_GET_STDDEV(&search_stats->unique_template_size);
}
uint64_t mapper_stats_template_length_get_expected_max(mapper_stats_t* const search_stats) {
  // Mean - 3 * StdDev
  const double mean = mapper_stats_template_length_get_mean(search_stats);
  const double max_dev = 6.0*mapper_stats_template_length_get_stddev(search_stats);
  return mean + max_dev;
}
uint64_t mapper_stats_template_length_get_expected_min(mapper_stats_t* const search_stats) {
  // Mean - 3 * StdDev
  const double mean = mapper_stats_template_length_get_mean(search_stats);
  const double max_dev = 6.0*mapper_stats_template_length_get_stddev(search_stats);
  return BOUNDED_SUBTRACTION(mean,max_dev,0.0);
}
double mapper_stats_template_length_get_sigma_dev(
    mapper_stats_t* const search_stats,
    const uint64_t template_length) {
  // Mean - 3 * StdDev
  const double mean = mapper_stats_template_length_get_mean(search_stats);
  const double dev = fabs((double)template_length-mean);
  const double stddev = mapper_stats_template_length_get_stddev(search_stats);
  return dev/stddev;
}

