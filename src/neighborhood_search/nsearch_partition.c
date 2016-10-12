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

#include "neighborhood_search/nsearch_partition.h"

/*
 * Compute Partition
 */
void nsearch_partition_compute(
    nsearch_partition_t* const nsearch_partition,
    const uint64_t chunk_offset,
    const uint64_t chunk_length) {
  // Pattern Partition
  nsearch_partition->offset_0 = chunk_offset;
  nsearch_partition->length_0 = chunk_length/2;
  nsearch_partition->offset_1 = chunk_offset + nsearch_partition->length_0;
  nsearch_partition->length_1 = chunk_length - nsearch_partition->length_0;
  // Single Region Chunk
  nsearch_partition->region_0 = NULL;
  nsearch_partition->region_1 = NULL;
}
void nsearch_partition_preconditioned_compute(
    nsearch_partition_t* const nsearch_partition,
    region_search_t* const regions,
    const uint64_t region_offset,
    const uint64_t num_regions) {
  // Region Partition
  nsearch_partition->region_offset_0 = region_offset;
  nsearch_partition->num_regions_0 = DIV_CEIL(num_regions,2);
  nsearch_partition->region_offset_1 = nsearch_partition->region_offset_0 + nsearch_partition->num_regions_0;
  nsearch_partition->num_regions_1 = num_regions - nsearch_partition->num_regions_0;
  // Search partition
  nsearch_partition->offset_0 = regions[nsearch_partition->region_offset_0].begin;
  const uint64_t last_region_idx_0 = nsearch_partition->region_offset_0+nsearch_partition->num_regions_0-1;
  nsearch_partition->length_0 = regions[last_region_idx_0].end - nsearch_partition->offset_0;
  nsearch_partition->offset_1 = regions[nsearch_partition->region_offset_1].begin;
  const uint64_t last_region_idx_1 = nsearch_partition->region_offset_1+nsearch_partition->num_regions_1-1;
  nsearch_partition->length_1 = regions[last_region_idx_1].end - nsearch_partition->offset_1;
  // Single Region Chunk
  nsearch_partition->region_0 = (nsearch_partition->num_regions_0==1) ? regions+nsearch_partition->region_offset_0 : NULL;
  nsearch_partition->region_1 = (nsearch_partition->num_regions_1==1) ? regions+nsearch_partition->region_offset_1 : NULL;
}
/*
 * Compute error partition
 */
void nsearch_partition_compute_error(
    nsearch_partition_t* const nsearch_partition,
    const uint64_t min_error,
    const uint64_t max_error) {
  // Parameters
  region_search_t* const region_0 = nsearch_partition->region_0;
  const uint64_t length_0 = nsearch_partition->length_0;
  region_search_t* const region_1 = nsearch_partition->region_1;
  const uint64_t length_1 = nsearch_partition->length_1;
  // Split error
  const uint64_t error_search = max_error/2;
  const uint64_t error_extend = max_error;
  /*
   * First partition (local)
   */
  // Base error
  nsearch_partition->search_0_min_error = 0;
  nsearch_partition->search_0_max_error = MIN(error_search,length_0); // Length constraint
  nsearch_partition->extend_1_local_min_error = 0;
  nsearch_partition->extend_1_local_max_error = MIN(error_extend,length_1); // Length constraint
  // Adjust using region-min/max
  if (region_0!=NULL) {
    nsearch_partition->search_0_min_error = region_0->min;
    nsearch_partition->search_0_max_error = MIN(nsearch_partition->search_0_max_error,region_0->max);
  }
  if (region_1!=NULL) {
    nsearch_partition->extend_1_local_min_error = region_1->min;
    nsearch_partition->extend_1_local_max_error = MIN(nsearch_partition->extend_1_local_max_error,region_1->max);
  }
  // Adjust using min-needed errors in each partition
  const uint64_t search_0_left_errors = max_error-nsearch_partition->extend_1_local_min_error;
  nsearch_partition->search_0_max_error = MIN(nsearch_partition->search_0_max_error,search_0_left_errors);
  const uint64_t extend_1_local_left_errors = max_error-nsearch_partition->search_0_min_error;
  nsearch_partition->extend_1_local_max_error = MIN(nsearch_partition->extend_1_local_max_error,extend_1_local_left_errors);
  /*
   * Second partition (local)
   */
  // Base error
  nsearch_partition->search_1_min_error = 0;
  nsearch_partition->search_1_max_error = MIN(error_search,length_1);
  nsearch_partition->extend_0_local_min_error = nsearch_partition->search_0_max_error+1;
  nsearch_partition->extend_0_local_max_error = MIN(error_extend,length_0);
  // Adjust using region min/max
  if (region_1!=NULL) {
    nsearch_partition->search_1_min_error = region_1->min;
    nsearch_partition->search_1_max_error = MIN(nsearch_partition->search_1_max_error,region_1->max);
  }
  if (region_0!=NULL) {
    nsearch_partition->extend_0_local_min_error = MAX(nsearch_partition->extend_0_local_min_error,region_0->min);
    nsearch_partition->extend_0_local_max_error = MIN(nsearch_partition->extend_0_local_max_error,region_0->max);
  }
  // Adjust combining constraints
  uint64_t search_1_left_errors = max_error-nsearch_partition->extend_0_local_min_error;
  nsearch_partition->search_1_max_error = MIN(nsearch_partition->search_1_max_error,search_1_left_errors);
  uint64_t extend_0_local_left_errors = max_error-nsearch_partition->search_1_min_error;
  nsearch_partition->extend_0_local_max_error = MIN(nsearch_partition->extend_0_local_max_error,extend_0_local_left_errors);
  /*
   * Region error (global)
   */
  if (region_0!=NULL) {
    nsearch_partition->extend_0_global_min_error = MAX(min_error,region_0->min);
    nsearch_partition->extend_0_global_max_error = max_error;
  } else {
    nsearch_partition->extend_0_global_min_error = min_error;
    nsearch_partition->extend_0_global_max_error = max_error;
  }
  if (region_1!=NULL) {
    nsearch_partition->extend_1_global_min_error = MAX(min_error,region_1->min);
    nsearch_partition->extend_1_global_max_error = max_error;
  } else {
    nsearch_partition->extend_1_global_min_error = min_error;
    nsearch_partition->extend_1_global_max_error = max_error;
  }
  // Adjust combining constraints
  nsearch_partition->extend_0_global_min_error =
      MAX(nsearch_partition->extend_0_global_min_error,nsearch_partition->extend_0_local_min_error);
  nsearch_partition->extend_1_global_min_error =
      MAX(nsearch_partition->extend_1_global_min_error,nsearch_partition->extend_1_local_min_error);
}
