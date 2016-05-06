/*
 * PROJECT: GEMMapper
 * FILE: neighborhood_search_hamming.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_PARTITION_H_
#define NSEARCH_PARTITION_H_

#include "utils/essentials.h"
#include "fm_index/fm_index.h"
#include "data_structures/interval_set.h"
#include "filtering/region_profile.h"

/*
 * Neighborhood-Search Error Partition
 */
typedef struct {
  // Pattern Regions
  uint64_t region_offset_0;
  uint64_t num_regions_0;
  uint64_t region_offset_1;
  uint64_t num_regions_1;
  // Single Region Chunk
  region_search_t* region_0;
  region_search_t* region_1;
  // Pattern Partition
  uint64_t offset_0;
  uint64_t length_0;
  uint64_t offset_1;
  uint64_t length_1;
  // Error Partition
  uint64_t search_0_min_error;
  uint64_t search_0_max_error;
  uint64_t search_1_min_error;
  uint64_t search_1_max_error;
  uint64_t extend_0_local_min_error;
  uint64_t extend_0_local_max_error;
  uint64_t extend_1_local_min_error;
  uint64_t extend_1_local_max_error;
  uint64_t extend_0_global_min_error;
  uint64_t extend_0_global_max_error;
  uint64_t extend_1_global_min_error;
  uint64_t extend_1_global_max_error;
} nsearch_partition_t;

/*
 * Compute Partition
 */
void nsearch_partition_compute(
    nsearch_partition_t* const nsearch_partition,
    const uint64_t chunk_offset,const uint64_t chunk_length);
void nsearch_partition_preconditioned_compute(
    nsearch_partition_t* const nsearch_partition,region_search_t* const regions,
    const uint64_t region_offset,const uint64_t num_regions);

/*
 * Compute error partition
 */
void nsearch_partition_compute_error(
    nsearch_partition_t* const nsearch_partition,
    const uint64_t min_error,const uint64_t max_error);

#endif /* NSEARCH_PARTITION_H_ */
