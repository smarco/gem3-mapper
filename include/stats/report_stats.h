/*
 * PROJECT: GEMMapper
 * FILE: report_stats.h
 * DATE: 06/06/2012
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

void init_mapping_stats(mapping_stats_t *mstats);

void collect_SE_mapping_stats(
    archive_search_t* const archive_search,
    matches_t* const matches,
    mapping_stats_t* mstats);
void collect_PE_mapping_stats(
    archive_search_t* const archive_search1,
    archive_search_t* const archive_search2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* mstats);

void output_mapping_stats(mapper_parameters_t *parameters,mapping_stats_t* mstats);
void merge_mapping_stats(mapping_stats_t *global_mstats,mapping_stats_t *mapping_stats,uint64_t num_threads);

#endif
