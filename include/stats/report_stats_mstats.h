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

typedef struct {
	 uint64_t reads[2][4];
	 uint64_t BSreads[2][2];
	 uint64_t unmapped[2];
	 uint64_t correct_pairs;
	 uint64_t base_counts[7][2][5];
	 uint64_t hist_mapq[256];
	 ihash_t *read_length_dist[2];
	 ihash_t *insert_size_dist;
	 ihash_t *distance_dist[2];
} mapping_stats_t;

#endif
