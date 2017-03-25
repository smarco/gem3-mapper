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

#define UNDERCONVERSION_CONTROL "NC_001416.1"
#define OVERCONVERSION_CONTROL "NC_001604.1"
#define SEQUENCING_CONTROL "NC_001422.1"

void init_mapping_stats(
    mapping_stats_t* const mstats);

void collect_se_mapping_stats(
    archive_search_t* const archive_search,
    matches_t* const matches,
    mapping_stats_t* mstats);
void collect_pe_mapping_stats(
    archive_search_t* const archive_search1,
    archive_search_t* const archive_search2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* mstats);

void output_mapping_stats(
    mapper_parameters_t* const parameters,
    mapping_stats_t* const mstats);
void merge_mapping_stats(
    mapping_stats_t* const global_mstats,
    mapping_stats_t* const mstats,
    const uint64_t num_threads);

#endif
