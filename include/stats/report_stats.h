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

#ifndef REPORT_STATS_H_
#define REPORT_STATS_H_

#include "system/commons.h"
#include "utils/hash.h"
#include "stats/report_stats_mstats.h"
#include "mapper/mapper.h"

#define UNDERCONVERSION_CONTROL "NC_001416.1:Lambda_5C_conversion"
#define OVERCONVERSION_CONTROL "NC_001604.1:T7_5mC_conversion"
#define SEQUENCING_CONTROL "NC_001422.1:PhiX"

void init_mapping_stats(
    mapping_stats_t *const mstats);

void collect_se_mapping_stats(
    const archive_search_t* archive_search,
    const matches_t* matches,
    mapping_stats_t* mstats);
void collect_pe_mapping_stats(
    const archive_search_t* archive_search1,
    const archive_search_t* archive_search2,
    paired_matches_t* paired_matches,
    mapping_stats_t* mstats);

void output_mapping_stats(
    const mapper_parameters_t* parameters,
    const mapping_stats_t* mstats);
void merge_mapping_stats(
    const mapping_stats_t* global_mstats,
    const mapping_stats_t* mstats,
    const uint64_t num_threads);

#endif
