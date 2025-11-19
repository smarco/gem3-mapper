/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Simon Heath  <simon.heath@gmail.com>
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
 * AUTHOR(S): Simon Heath  <simon.heath@gmail.com>
 */

#ifndef REPORT_STATS_MSTATS_H_
#define REPORT_STATS_MSTATS_H_

#include "system/commons.h"
#include "utils/hash.h"

#define N_BASE_COUNTS 5

/*
 * Index 1: read end
 * Index 2: base (N, A, C, G, T) 
 */
typedef struct {
    uint64_t counts[2][N_BASE_COUNTS];
} base_counts_t;

void clear_base_counts(base_counts_t *p);
void add_base_counts(base_counts_t *p, base_counts_t const *p1);

typedef struct {
     int n_control_seq;
	 uint64_t *reads;
	 uint64_t BSreads[2][2];
	 uint64_t unmapped[2];
	 uint64_t correct_pairs;
	
	 /* Base counts */
	 base_counts_t overall_counts;
	 /*
	  * Base counts for different categories
	  * Index: C2T, G2A
	  *
	  * First entries are for general, followed by the conversion controls  
	  */
	 base_counts_t* base_counts[2];

	 // uint64_t base_counts[7][2][5];
	 uint64_t hist_mapq[256];
	 ihash_t *read_length_dist[2];
	 ihash_t *insert_size_dist;
	 ihash_t *distance_dist[2];
} mapping_stats_t;

#define mapping_stats_reads_get(mstats, read, rix) (mstats->reads + (rix * 2 + read))
#define mapping_stats_reads_add(mstats, read, rix, x) ((*mapping_stats_reads_get(mstats, read, rix)) += x)
#define mapping_stats_reads_inc(mstats, read, rix) ((*mapping_stats_reads_get(mstats, read, rix))++)

#define mapping_stats_base_counts(mstats, bs_strand, idx) (mstats->base_counts[bs_strand] + idx)

mapping_stats_t *new_mapping_stats(int n_control_seq);
void setup_mapping_stats(mapping_stats_t *ms, int n_control_seq);

#endif
