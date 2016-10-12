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
 */

#ifndef PAIRED_MATCHES_H_
#define PAIRED_MATCHES_H_

#include "utils/essentials.h"
#include "archive/search/archive_search_pe_parameters.h"
#include "text/sequence.h"
#include "matches/matches.h"
#include "matches/paired_map.h"
#include "mapper/mapper_stats.h"

/*
 * Paired-Matches Classes
 */
typedef enum {
  paired_matches_class_unmapped = 0,
  paired_matches_class_tie_perfect = 1,
  paired_matches_class_tie = 2,
  paired_matches_class_mmap_d1 = 3,
  paired_matches_class_mmap = 4,
  paired_matches_class_unique = 5,
} paired_matches_class_t;
extern const char* paired_matches_class_label[6];

/*
 * Paired Matches
 */
typedef struct {
  /* State */
  paired_matches_class_t paired_matches_class; // Matches Class
  uint64_t max_complete_stratum;               // MCS
  /* Matches Counters */
  matches_counters_t* counters;                // Counters
  /* Single-End Matches */
  matches_t* matches_end1;                     // Matches end1
  matches_t* matches_end2;                     // Matches end2
  /* Paired-End Matches */
  vector_t* paired_maps;                       // Paired Maps (paired_map_t)
  vector_t* discordant_paired_maps;            // Discordant Paired Maps (paired_map_t)
  /* Paired-Matches Metrics */
  matches_metrics_t metrics;                   // PE-Matches Metrics
} paired_matches_t;

/*
 * Setup
 */
paired_matches_t* paired_matches_new(void);
void paired_matches_clear(paired_matches_t* const paired_matches,const bool clear_matches);
void paired_matches_delete(paired_matches_t* const paired_matches);

/*
 * Accessors
 */
bool paired_matches_is_mapped(const paired_matches_t* const paired_matches);
uint64_t paired_matches_get_num_maps(const paired_matches_t* const paired_matches);
paired_map_t* paired_matches_get_maps(paired_matches_t* const paired_matches);

uint64_t paired_matches_counters_get_count(paired_matches_t* const paired_matches,const uint64_t distance);
uint64_t paired_matches_counters_get_total_count(paired_matches_t* const paired_matches);

uint64_t paired_matches_get_first_stratum_matches(paired_matches_t* const paired_matches);

/*
 * Adding Paired-Matches
 */
void paired_matches_add(
    paired_matches_t* const paired_matches,
    match_trace_t* const match_trace_end1,
    match_trace_t* const match_trace_end2,
    const pair_relation_t pair_relation,
    const pair_orientation_t pair_orientation,
    const pair_layout_t pair_layout,
    const uint64_t template_length,
    const double template_length_sigma);

/*
 * Finding Pairs
 */
void paired_matches_find_pairs(
    paired_matches_t* const paired_matches,
    const search_paired_parameters_t* const search_paired_parameters,
    mapper_stats_t* const mapper_stats);
void paired_matches_find_discordant_pairs(
    paired_matches_t* const paired_matches,
    const search_paired_parameters_t* const search_paired_parameters);

/*
 * Filters
 */
void paired_matches_filter_by_mapq(
    paired_matches_t* const paired_matches,
    const uint8_t mapq_threshold);

/*
 * Sort
 */
void paired_matches_sort_by_swg_score(paired_matches_t* const paired_matches);

/*
 * Display
 */
void paired_matches_print(FILE* const stream,paired_matches_t* const paired_matches);

#endif /* PAIRED_MATCHES_H_ */
