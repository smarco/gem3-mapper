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
} mapping_stats_t;

#endif
